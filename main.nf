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
include { infereIBD } from './modules/snipar.nf'
include { imputeGenotypes } from './modules/snipar.nf'
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
    /* INPUT FILES */
    segmentByChromosome(
        params.bed,
        params.bim,
        params.fam
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

    /* SNIPAR 
        The .out.map { i -> i[0].parent } processes are necessary
        because SNIPAR looks for all files within a directory.
        See: https://snipar.readthedocs.io/en/latest/tutorial.html#tutorial
    */
    infereIBD(
        segmentByChromosome.out.map { i -> i[0].parent },
        createPedigree.output.pedigree,
        params.chr_range
    )
    imputeGenotypes(
        infereIBD.out.map { i -> i[0].parent },
        segmentByChromosome.out.map { i -> i[0].parent },
        createPedigree.output.pedigree,
        params.chr_range
    )
    runSiblingsFGWAS(
        createPhenotype.output.phenotype,
        segmentByChromosome.out.map { i -> i[0].parent },
        imputeGenotypes.out.map { i -> i[0].parent },
        createPedigree.output.pedigree,
        createAgesexKinship.output.kinship,
        params.chr_range
    )
}