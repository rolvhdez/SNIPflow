#!/usr/bin/env nextflow

process runSiblingsFGWAS {
    publishDir "${params.outDir}/sumstats", mode: 'copy'

    input:
        tuple path(bed), path(bim), path(fam)
        path phenotype
        path pedigree
        val chr

    output:
        path "chr_${chr}.sumstats.gz", emit: "sumstats"
        path "chr_${chr}.hdf5", emit: "hdf5"
        path "fgwas_${chr}.log", emit: "log"

    script:
    def bedfile = bed.baseName
    """
    #!/bin/bash

    gwas.py ${phenotype} \
        --bed ${bedfile} \
        --pedigree ${pedigree} \
        --cpu 4 --threads 1 \
        --out chr_${chr}
    """
}