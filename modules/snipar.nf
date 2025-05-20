#!/usr/bin/env nextflow

process infereIBD {
    containerOptions '--user root'
    input:
        path segmentedDir
        path pedigree
        val chr_range

    output:
        path "chr_*.ibd.segments.gz", emit: "ibd_rel"

    script:
    def bedfile = "${segmentedDir}/chr_@"
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

process imputeGenotypes {
    publishDir "${params.outDir}/impute", mode: 'copy'
    containerOptions '--user root'
    input:
    path segmentedDir
    path ibdDir
    path pedigree
    val chr_range

    output:
    path "chr_*.imputed.hdf5", emit: "chr_imputed"

    script:
    def bedfile = "${segmentedDir}/chr_@"
    def ibd = "${ibdDir}/chr_@.ibd"

    """
    #!/bin/bash

    impute.py \
        --ibd "${ibd}" \
        --bed "${bedfile}" \
        --pedigree "${pedigree}" \
        --chr_range "${chr_range}" \
        --threads 10 \
        --out chr_@.imputed \
    """
}

process runSiblingsFGWAS {
    publishDir "${params.outDir}/sumstats", mode: 'copy'
    containerOptions '--user root'

    input:
    path segmentedDir
    path imputedDir
    path phenotype
    path pedigree
    path kinship
    val chr_range

    output:
    path "chr_*_robust.sumstats.gz", emit: "sumstats"
    path "chr_*_robust.hdf5", emit: "hdf5"

    script:
    def bedfile = "${segmentedDir}/chr_@"
    def imputed = "${imputedDir}/chr_@.imputed"
    def ibdrel = "${kinship.baseName}"
    """
    #!/bin/bash

    gwas.py ${phenotype} \
        --phen_index 1 \
        --bed ${bedfile} \
        --imp ${imputed} \
        --pedigree ${pedigree} \
        --chr_range ${chr_range} \
        --ibdrel_path ${ibdrel} --sparse_thresh 0.1 \
        --robust \
        --cpus 10 --threads 4 \
        --out chr_@_robust
    """
}