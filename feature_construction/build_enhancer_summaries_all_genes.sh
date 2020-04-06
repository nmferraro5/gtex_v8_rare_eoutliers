#!/bin/bash

# author: Emily Tsang

# compile the site-level feature for each genes
# only include features for individuals that have at least one variant within the vicinity of the gene

set -o nounset -o errexit -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

# takes as input several parameters:
# * base directory in which to write; will create this directory with subdirectories
# * window size
# * bed file with intervals around which to build the window (either with TSSs or entire gene bodies)
# * bed file with coding regions to exclude (optional)

### SET THE NUMBER OF PROCESSES TO RUN IN PARALLEL
nproc=10

if [ $# -ge 4 ] && [ $# -le 5 ]; then
    indir=$1
    outdir=$2
    windowsize=$3
    gene=$4
    if [ $# -eq 5 ]; then
	coding=$5
    else
	# make dummy file that is empty
	if [ ! -e dummy.bed ]; then
	    echo -e "chr1\t0\t1" > dummy.bed
	fi
	coding=dummy.bed
    fi
else
    echo "usage: build_enhancer_summaries_all_genes.sh v7/multi_omics input_directory output_directory window_size gene_bed (exclude_bed)"
    exit
fi

# make output directory (and subdirectories) if not already there
# if the upper-level directory is there, assume the subdirectories are also there
current_dir=`pwd` # manually added a MAF0-double file
if [ ! -d $outdir ]; then
    mkdir -p $outdir
    cd $outdir
    mkdir -p MAF0-1
    mkdir -p MAF1-5
    mkdir -p MAF5-10
    mkdir -p MAF10-25
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

# get window of variants around each gene and feed that into feature building script
# do this for different MAF thresholds
# MAF 0-1%, 1-5%, 5-10%, 10-25% where the beginning of the range is excluded and the end is included
makegenefeatures() {
    var_features=$1
    filter1kg=$2

    #endpoints=(0 0.01 0.05 0.1 0.25)
    #names=(0 1 5 10 25)
    endpoints=(0 0.01)
    names=(0 1)
    #endpoints=(0 0.00141) #0.00141 for SNPs/indels and 0.0019 for SV
    #endpoints=(0 0.00071) #0.00141 for SNPs/indels and 0.0019 for SV
    #names=(0 single)

    # strip features bed file name to get individual id and variant type
    fname=`basename $var_features`
    sample=${fname%%_*}
    vartype=${fname%_promoters.bed.gz}
    vartype=${vartype#*_}
    # create output file name
    out_prefix=${sample}_promoters_bygene
    
    # compile features at the gene level
    # can check if the MAF is also rare in the 1KG super population 
    pop=0
    for ((i=1; i<${#endpoints[@]}; i++))
    do
	mafname=${names[$i-1]}-${names[$i]}
	zcat $var_features | cut -f1-18 | \
            bedtools sort | \
            awk -v start=${endpoints[$i-1]} -v end=${endpoints[$i]} -v filter=$filter1kg -v offset=${pop} '$4>start && $4<=end && (filter!="T" || end>0.01 || $(NF-offset)=="NA" || $(NF-offset) <=0.001)' | \
	    bedtools intersect -sorted -wa -v -a stdin -b $coding | \
	    bedtools window -r $windowsize -l $windowsize -a $gene -b stdin | \
	    python ${scriptdir}/build_enhancer_summaries.py --features ${outdir}/MAF${mafname}/${vartype}/${out_prefix}.txt
            #cut -f4-38 > ${outdir}/MAF${mafname}/${vartype}/${out_prefix}.txt
    done
}

# make things available to parallel (note: cannot export arrays, so endpoints and names are defined in the function)
export -f makegenefeatures
export scriptdir
export coding
export windowsize
export gene
export outdir

parallel --jobs $nproc makegenefeatures ::: ${indir}/*_enhancers.bed.gz ::: T

echo "all done!"

