#!/bin/bash

# Author: Emily Tsang
# Updated: Nicole Ferraro, summer 2018

set -o nounset -o errexit -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

## run add_features_variant_beds.sh in batches

## arguments: takes as input the version and (optionally) the number of cores to use
if [ $# -eq 0 ]; then
    ncores=25
elif [ $# -eq 1 ]; then
    ncores=$1
else
    echo "usage: run_add_features_variant_beds.sh [ncores]"
    exit
fi

## set output directory
dir=${RAREDIR}/features_v8

outdir=${dir}/bySite

mkdir -p $outdir

parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individualsRA/*_SNPs.bed.gz ::: ${KG_AF_SNPS} ::: ${GNOMAD_SNPS} ::: ${outdir}

parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_indels.bed.gz ::: ${KG_AF_INDELS} ::: ${GNOMAD_INDELS} ::: ${outdir}

parallel --jobs $ncores --xapply ${scriptdir}/add_features_variant_beds.sh ::: ${dir}/variantBeds/individuals/*_HallLabSV.bed.gz ::: ${AF_SVs} ::: ${outdir}

wait
