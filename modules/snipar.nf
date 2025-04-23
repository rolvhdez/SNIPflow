#!/usr/bin/env nextflow

process infereIBD {
    publishDir "${params.outDir}/ibd_inference", mode: 'copy'
    containerOptions '--user root'
    input:
    path bed_path
    path pedigree
    val chr_range

    output:
    path "chr_*.ibd.segments.gz", emit: "ibd_rel"

    script:
    """
    #!/bin/bash

    ibd.py \
        --bed ${bed_path}/chr_@ \
        --pedigree ${pedigree} \
        --chr_range ${chr_range} \
        --threads 4 --ld_out --batches 4 \
        --out chr_@
    """
}

process imputeGenotypes {
    input:
    path ibd_path
    path bed_path
    path pedigree
    val chr

    output:
    path "chr_*.imputation.hdf5", emit: "imputed"

    script:
    """
    #!/bin/bash

    impute.py \
        --ibd ${ibd_path}/chr_@.ibd \
        --bed ${bed_path}/chr_@ \
        --pedigree ${pedigree} \
        --chr_range ${chr} \
        --threads 10 \
        --out chr_@.imputation
    """
}

process runSiblingsFGWAS {
    publishDir "${params.outDir}/sumstats", mode: 'copy'
    containerOptions '--user root'

    input:
    path phenotype
    path bed_path
    path imp_path
    path pedigree
    path kinship
    val chr_range

    output:
    path "chr_*.sumstats.gz", emit: "sumstats"
    path "chr_*.hdf5", emit: "hdf5"

    script:
    """
    #!/bin/bash

    gwas.py ${phenotype} \
        --bed ${bed_path}/chr_@ \
        --imp ${imp_path}/chr_@ \
        --pedigree ${pedigree} \
        --ibdrel_path ${kinship.baseName} --sparse_thresh 0.1 \
        --phen_index 1 \
        --robust \
        --cpus 10 --threads 10 \
        --chr_range ${chr_range} \
        --out chr_@
    """
}