#!/usr/bin/env nextflow

process infereIBD {
    containerOptions '--user root'
    input:
        path segmentedFiles
        path pedigree
        val chr_range

    output:
        path "chr_*.ibd.segments.gz", emit: "ibd_rel"

    script:
    def bedfile = "${segmentedFiles}/chr_@"
    """
    #!/bin/bash

    ibd.py \
        --bed ${bedfile} \
        --pedigree ${pedigree} \
        --chr_range ${chr_range} \
        --threads 4 --ld_out --batches 4 \
        --out chr_@
    """
}

process runSiblingsFGWAS {
    publishDir "${params.outDir}/sumstats", mode: 'copy'
    containerOptions '--user root'

    input:
//        tuple path(bed), path(bim), path(fam)
        path segmentedFiles
        path phenotype
        path pedigree
        path kinship
        val chr

    output:
        path "chr_*.sumstats.gz", emit: "sumstats"
        path "chr_*.hdf5", emit: "hdf5"

    script:
    def bedfile = "${segmentedFiles}/chr_@"
    def ibdrel = "${kinship.baseName}"
    """
    #!/bin/bash

    gwas.py ${phenotype} \
        --phen_index 1 \
        --bed ${bedfile} \
        --pedigree ${pedigree} \
        --chr_range ${chr} \
        --ibdrel_path ${ibdrel} --sparse_thresh 0.05 \
        --cpu 8 --threads 1 \
        --out chr_@
    """
}