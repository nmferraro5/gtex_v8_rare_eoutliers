#!/bin/bash

# Author: Emily Tsang
# Updated: Nicole Ferraro, summer 2018

set -o nounset -o errexit -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

## run add_enhancers_variant_beds.sh in batches

## arguments: takes as input the version and (optionally) the number of cores to use
if [ $# -eq 0 ]; then
    ncores=10
elif [ $# -eq 1 ]; then
    ncores=$2
else
    echo "usage: run_add_enhancers_variant_beds.sh [ncores]"
    exit
fi

## set output directory
dir=${RAREDIR}/features_v8

outdir=${dir}/bySite
AF_SVS=${AF_SVs}
mkdir -p $outdir

parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_SNPs.bed.gz ::: ${GNOMAD_SNPS} ::: ${outdir}

parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_indels.bed.gz ::: ${GNOMAD_INDELS} ::: ${outdir}

outdir=${dir}/bySite/HallLabSV
parallel --jobs $ncores --xapply ${scriptdir}/add_predEnhancers_variant_beds.sh ::: ${dir}/variantBeds/individuals/HallLabSV_hg38/*_HallLabSV.bed.gz ::: ${AF_SVS} ::: ${outdir}

wait
