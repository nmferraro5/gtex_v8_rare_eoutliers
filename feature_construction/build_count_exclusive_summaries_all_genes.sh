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
    #indir_SV=$3
    outdir=$2
    windowstart=$3
    windowend=$4
    gene=$5
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
if [ -d ${outdir}/counts ]; then # remove ! for now to create new directories
    cd $outdir
    mkdir -p counts
    cd counts
    mkdir -p MAF0-double
    mkdir -p MAF0-single
    mkdir -p MAFnovel
    mkdir -p MAFgnomad0-1only
    mkdir -p MAFgnomad0-1
    mkdir -p MAFgnomad1-5
    for dir in MAF*
    do
	cd $dir
	mkdir -p SNPs
	mkdir -p indels
	mkdir -p HallLabSV
	cd ..
    done
fi
cd $current_dir

# Start separately if including genebody (windowstart = 0)
if [ $windowstart -eq 0 ]; then
    cat $gene | \
        awk -v lb=$windowend 'BEGIN {OFS="\t"} { print $1,$2-lb,$3+lb,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($3>0) print $1,$2,$3,$4; else print $1,$2,0,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($2!=0 && $3!=0) print $1,$2,$3,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($1 != "chrM") print $1,$2,$3,$4}' > genetemp.txt
else
    bedtools flank -i $gene -g ${RAREDIR}/preprocessing_v8/chrom_length_dummy.txt -b $windowend | \
        awk -v lb=$windowstart 'BEGIN {OFS="\t"} { if (NR % 2 == 1) print $1,$2,$3-lb,$4; else print $1,$2+lb,$3,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($3>0) print $1,$2,$3,$4; else print $1,$2,0,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($2!=0 && $3!=0) print $1,$2,$3,$4}' | \
        awk 'BEGIN {OFS="\t"} { if ($1 != "chrM") print $1,$2,$3,$4}' > genetemp.txt
fi

# get window of variants around each gene and count them
# do this for different MAF thresholds
# MAF 0-1%, 1-5%, 5-10%, 10-25% where the beginning of the range is excluded and the end is included

# takes as input a variant bed file and whether or not to apply the gnomad filter
# the second argument must be T to apply gnomad filter (otherwise the filter will be applied)
# if the gnomad filter is to be applied, the second to last column of the file must have the gnomad EUR AF
makecountfeatures() {
    var_bed=$1
    filter1kg=$2

    # strip features bed file name to get gtex id and variant type
    echo "var_bed is " $var_bed
    fname=`basename $var_bed`
    echo "fname is " $fname
    sample=${fname%%_*}
    echo "sample is " $sample
    vartype=${fname%.bed.gz}
    vartype=${vartype#*_}
    vartype=${vartype%_*}
    
    out_prefix=${sample}_counts_bygene

    # compile features at the gene level
    # can check if the MAF is also rare in the 1KG super population 
    pop=0

    ## Just for HallLabSV
    mafname=MAF0-single
    outfile1=${outdir}/counts/MAF0-single/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile1
    zcat $var_bed | \
        awk -v offset=${pop} '$4<=0.00095' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile1
    wait
    mafname=MAF0-1
    outfile2=${outdir}/counts/MAF0-1/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile2
    zcat $var_bed | \
        awk -v offset=${pop} '$4<=0.01' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile2
    wait
    mafname=MAF1-5
    outfile2=${outdir}/counts/MAF1-5/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile2
    zcat $var_bed | \
        awk -v offset=${pop} '$4>0.01 && $4<=0.05' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile2
    mafname=novel
    outfile1=${outdir}/counts/MAFnovel/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile1
    zcat $var_bed | \
        awk -v offset=${pop} '$(NF-offset)=="NA" || $(NF-offset) == 0' | \
        awk -v offset=${pop} '$(NF-offset) !~ /,/ { print }' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile1
    wait
    mafname=gnomad0-1only
    outfile4=${outdir}/counts/MAFgnomad0-1only/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile4
    zcat $var_bed | \
        awk -v offset=${pop} '$(NF-offset) <= 0.01' | \
        awk -v offset=${pop} '$(NF-offset) !~ /,/ { print }' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile4
    zcat $var_bed | \
        awk -v offset=${pop} '$(NF-offset)=="NA" || $(NF-offset) <= 0.01' | \
        awk -v offset=${pop} '$(NF-offset) !~ /,/ { print }' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile4
    wait
    mafname=gnomad1-5
    outfile5=${outdir}/counts/MAFgnomad1-5/${vartype}/${out_prefix}.txt
    echo -e "gene_id\tn_variants" > $outfile5
    zcat $var_bed | \
        awk -v offset=${pop} '$(NF-offset) > 0.01 && $(NF-offset)>= 0.05' | \
        awk -v offset=${pop} '$(NF-offset) !~ /,/ { print }' | \
        bedtools intersect -sorted -wa -v -a stdin -b $coding | \
        bedtools intersect -c -a genetemp.txt -b stdin | \
        bedtools groupby -i stdin -g 4 -c 5 >> $outfile5
   
    ## create output file name
    out_prefix=${sample}_counts_bygene
    
    # compile features at the gene level
    pop=0

    for ((i=1; i<${#endpoints[@]}; i++))
    do
        mafname=${names[$i-1]}-${names[$i]}
        outfile=${outdir}/counts/MAF${mafname}/${vartype}/${out_prefix}.txt
        echo -e "gene_id\tn_variants" > $outfile
        zcat $var_bed | \
            awk -v start=${endpoints[$i-1]} -v end=${endpoints[$i]} -v filter=$filter1kg -v offset=${pop} '$4>start && $4<=end && (filter!="T" || end>0.01 || $(NF-offset)=="NA" || $(NF-offset) <=0.01)' | \
            bedtools intersect -sorted -wa -v -a stdin -b $coding | \
            bedtools intersect -c -a genetemp.txt -b stdin | \
            bedtools groupby -i stdin -g 4 -c 5 >> $outfile
    done
    echo $outfile
}


# make things available to parallel (note: cannot export arrays, so endpoints and names are defined in the function)
export -f makecountfeatures
export coding
export windowstart
export windowend
export gene
export outdir

echo "Processing SNPs and indels..."
parallel --jobs $nproc makecountfeatures ::: ${indir}/*_features.bed.gz ::: T

echo "Processing SVs..."
svdir=${RAREDIR}/features_v7/variantBeds/individuals/HallLabSV_hg38
parallel --jobs $nproc makecountfeatures ::: ${svdir}/GTEX-*_HallLabSV.bed.gz ::: F

echo "all done!"
