#!/usr/bin/env bash
docker run \
    --rm \
    -v /shared/:/shared/ \
    --memory="10G" \
    --cpus="2" \
    -t carynwillis94/coloc-r:latest \
    Rscript /shared/shared/cdwillis/cancer_GWAS/coloc_analysis/coloc_step.R $1 $2 $3

