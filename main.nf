#!/usr/bin/env nextflow

params.baseline = null
params.kinship = null
params.bed = null
params.bim = null
params.fam = null
params.chr_range = "1-22"
params.outDir = "./results"

/* Processes */
include { createAgesexKinship } from './modules/input-preps.nf'
include { createPedigree } from './modules/input-preps.nf'
include { createPhenotype } from './modules/input-preps.nf'
include { segmentByChromosome } from './modules/input-preps.nf'
include { runSiblingsFGWAS } from './modules/snipar.nf'

/* Functions */	
def expandRanges(String str) {
    // Split the CHR_RANGE string into individual elements
    def elements = str.split()
    def chrom = []
        elements.each { element ->
            if (element.contains('-')) {
                def (start, end) = element.split('-').collect { it as int }
                chrom.addAll(start..end)
            } else {
                chrom.add(element as int)
            }
        }
    return chrom
}

workflow {
    // Find the grouped PLINK files (.bed, .bim, .fam)
//    Channel
//        .fromFilePairs ("${params.genotype}/*.{bed,fam,bim}", size:3, flat: true)
//        .ifEmpty { error "No matching plink files" }
//        .set { raw_plink_data }

    // Make a list of the chromosomes to use
    Channel
        .of(expandRanges(params.chr_range.toString()))
        .ifEmpty { error "No chromosomes found in the range ${params.chr_range}" }
        .flatten()
        .set { chr_channel }

    /* INPUT FILES */
    segmentByChromosome(
        params.bed,
        params.bim,
        params.fam,
        chr_channel
    )
    createAgesexKinship(
        params.baseline,
        params.kinship
    )
    createPedigree(
        createAgesexKinship.output.agesex, 
        createAgesexKinship.output.kinship
    )
    createPhenotype(
        params.baseline,
        params.kinship,
        createPedigree.output
    )

    /* SNIPAR */
    runSiblingsFGWAS(
        segmentByChromosome.output,
        createPhenotype.output,
        createPedigree.output,
        chr_channel
    )
}