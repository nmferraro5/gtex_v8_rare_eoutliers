#!/bin/bash

# author: Emily Tsang

# count variants within vicinity of each gene for every individual

set -o nounset -o errexit -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

# takes as input several parameters:
# * base directory in which to write; will create this directory with subdirectories
# * window size
# * bed file with intervals around which to build the window (either with TSSs or entire gene bodies)
# * bed file with coding regions to exclude (optional)

# note that rare variants are required to also be rare in the 1KG EUR population

### SET THE NUMBER OF PROCESSES TO RUN IN PARALLEL
nproc=15

if [ $# -ge 4 ] && [ $# -le 6 ]; then
    indir=$1
    #indir_SV=$2
    outdir=$2
    windowsize=$3
    gene=$4
    if [ $# -eq 6 ]; then
	coding=$6
    else
	# make dummy file that is empty
	if [ ! -e dummy.bed ]; then
	    echo -e "chr1\t0\t1" > dummy.bed
	fi
	coding=dummy.bed
    fi
else
    echo "usage: build_count_summaries_all_genes.sh site_features_bed output_directory window_size gene_bed (exclude_bed)"
    exit
fi

# make output directory (and subdirectories) if not already there
# if the counts subdirectory is there, assume the nested subdirectories are also there
current_dir=`pwd`
mkdir -p $outdir
if [ ! -d ${outdir}/counts ]; then
    cd $outdir
    mkdir counts
    cd counts
    mkdir MAF0-1
    mkdir MAF1-5
    mkdir MAF5-10
    mkdir MAF10-25
    for dir in MAF*
    do
	cd $dir
	mkdir SNPs
	mkdir indels
	mkdir HallLabSV
	cd ..
    done
fi
cd $current_dir
# get window of variants around each gene and count them
# do this for different MAF thresholds
# MAF 0-1%, 1-5%, 5-10%, 10-25% where the beginning of the range is excluded and the end is included

# takes as input a variant bed file and whether or not to apply the 1KG filter
# the second argument must be T to apply 1KG filter (otherwise the filter will be applied)
# if the 1KG filter is to be applied, the second to last column of the file must have the 1KG EUR MAFs
makecountfeatures() {
    var_bed=$1
    filter1kg=$2

    endpoints=(0 0.01 0.05 0.1 0.25)
    names=(0 1 5 10 25)
    #endpoints=(0 0.01)
    #names=(0 1)
    #endpoints=(0 0.001 0.01)
    #names=(0 0.1 1)

    # strip features bed file name to get gtex id and variant type
    echo "var_bed is " $var_bed
    fname=`basename $var_bed`
    echo "fname is " $fname
    sample=${fname%%_*}
    echo "sample is " $sample
    vartype=${fname%.bed.gz}
    vartype=${vartype#*_}
    vartype=${vartype%_*}
   
    # create output file name
    out_prefix=${sample}_counts_bygene
    
    # compile features at the gene level
    # can check if the MAF is also rare in the 1KG super population 
    pop=0
    for ((i=1; i<${#endpoints[@]}; i++))
    do
	mafname=${names[$i-1]}-${names[$i]}
	outfile=${outdir}/counts/MAF${mafname}/${vartype}/${out_prefix}.txt
	echo -e "gene_id\tn_variants" > $outfile
	zcat $var_bed | \
	    bedtools sort | \
	    awk -v start=${endpoints[$i-1]} -v end=${endpoints[$i]} -v filter=$filter1kg -v offset=${pop} '$4>start && $4<=end && (filter!="T" || end>0.01 || $(NF-offset)=="NA" || $(NF-offset) <=0.01)' | \
	    bedtools intersect -sorted -wa -v -a stdin -b $coding | \
	    bedtools window -c -r $windowsize -l $windowsize -a $gene -b stdin | \
	    cut -f4,5 \
	    >> $outfile
    done
    echo $outfile
}

makevariantlist() {
    var_bed=$1
    filter1kg=$2
    outfile=$3

    fname=`basename $var_bed`
    sample=${fname%%_*}

    pop=0
    echo $var_bed
    zcat $var_bed | \
	bedtools sort | \
        awk -v start=0 -v end=0.01 -v filter=$filter1kg -v offset=${pop} '$4>start && $4<=end && (filter!="T" || $(NF-offset)=="NA" || $(NF-offset) <=0.01)' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools window -r $windowsize -l $windowsize -a $gene -b stdin | \
        awk -v ind=$sample 'BEGIN{OFS="\t"}{print ind,$4,$5,$6,$7,$8}' \
        >> $outfile
}

# make things available to parallel (note: cannot export arrays, so endpoints and names are defined in the function)
export -f makecountfeatures
export coding
export windowsize
export gene
export outdir

echo "Processing SNPs and indels..."
parallel --jobs $nproc makecountfeatures ::: ${indir}/*_features.bed.gz ::: T

echo "Processing SVs..."
svdir=${RAREDIR}/features_v7/variantBeds/individuals/HallLabSV_hg38
parallel --jobs $nproc makecountfeatures ::: ${svdir}/GTEX-*_HallLabSV.bed.gz ::: F

echo "Processing all rare variants together into a single file for each variant type..."
# do NOT run this in parallel because that doesn't play nice with appending to the same file
# first erase contents of all_rare_out files (or create them if they don't exist)
all_rare_out=${outdir}/all_rare_variants
all_out_files="${all_rare_out}_SNPs.txt ${all_rare_out}_indels.txt ${all_rare_out}_HallLabSV.txt"
for f in $all_out_files
do
    echo -n "" > $f
done

echo "SNPs..."
##loop through files one by one
snpout=${all_rare_out}_SNPs.txt
for snp in ${indir}/*_SNPs_features.bed.gz
do
    makevariantlist $snp T $snpout
done

echo "indels..."
indelout=${all_rare_out}_indels.txt
for indel in ${indir}/*_indels_features.bed.gz
do
    makevariantlist $indel T $indelout
done

if [ $svdir != "NONE" ]; then
    echo "SVs..."
    svout=${all_rare_out}_HallLabSV.txt
    for sv in  ${svdir}/GTEX-*_HallLabSV.bed.gz
    do
	makevariantlist $sv F $svout
    done
else
    rm ${all_rare_out}_HallLabSV.txt
fi
echo "all done!"
