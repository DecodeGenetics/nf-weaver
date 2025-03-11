#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/**
 * Nextflow modules
 */
include { WEAVER_INDEX } from  './modules/weaver_index.nf'\
    addParams(WEAVER_INDEX_PUBLISH_DIR: "${launchDir}/weaver_index",\
              WEAVER_INDEX_PUBLISH_MODE: "symlink",\
              WEAVER_INDEX_LABEL: "weaver_index")

include { WEAVER_MAP } from './modules/weaver_map.nf'\
    addParams(WEAVER_MAP_PUBLISH_DIR: "${launchDir}/weaver_map",\
              WEAVER_MAP_PUBLISH_MODE: "symlink",\
              WEAVER_MAP_LABEL: "weaver_map")

include { WEAVER_REMAP } from './modules/weaver_remap.nf'\
    addParams(WEAVER_REMAP_PUBLISH_DIR: "${launchDir}/weaver_remap",\
              WEAVER_REMAP_PUBLISH_MODE: "symlink",\
              WEAVER_REMAP_LABEL: "weaver_remap")

include { DEEP_VARIANT } from './modules/deep_variant.nf'
include { GRAPHTYPER   } from './modules/graphtyper.nf'

/**
 * Parameter validation
 */
if (params.graph == "") {
    println "ERROR: No provided '--graph <gfa>' to map to. See README for details."
    System.exit(1)
}

if (params.genome == "") {
    println "ERROR: No '--genome <fa>' specified. See README for details."
    System.exit(1)
}

if (params.map == "" && params.remap == "") {
    println "ERROR: Neither '--map <tsv>' nor '--remap <tsv>' specified. See README for details."
    System.exit(1)
}

/**
 * Workflow variables
 */
def graph = file(params.graph, checkIfExists: true)
def graph_wmi = file(params.graph + ".wmi", checkIfExists: false)
def make_index = !graph_wmi.exists()
def indexed_genome = [params.genome, params.genome + ".fai"]
def indexed_graph = [graph, graph_wmi]

/**
 * Weaver remap workflow
 */
workflow {
    if (make_index) {
        WEAVER_INDEX(graph, params.weaver_k, params.weaver_w)
        indexed_graph = WEAVER_INDEX.out.indexed_graph
    }

    // Map jobs
    if (params.map != "") {
        Channel.fromPath(params.map)
            .splitCsv(sep:"\t", header:false)
            .set { bams_to_map }

        WEAVER_MAP(indexed_graph, params.genome, bams_to_map)

        // Deep variant
        if (params.run_deep_variant) {
            DEEP_VARIANT(WEAVER_MAP.out.sample_name, indexed_genome, WEAVER_MAP.out.indexed_cram)
        }

        // Graphtyper
        if (params.run_graphtyper) {
            GRAPHTYPER(WEAVER_MAP.out.sample_name, indexed_genome, WEAVER_MAP.out.indexed_cram)
        }
    }

    // Remap jobs
    if (params.remap != "") {
        Channel.fromPath(params.remap)
            .splitCsv(sep:"\t", header:false)
            .set { bams_to_remap }

        WEAVER_REMAP(indexed_graph, params.genome, bams_to_remap)

        // Deep variant
        if (params.run_deep_variant) {
            DEEP_VARIANT(WEAVER_REMAP.out.sample_name, indexed_genome, WEAVER_REMAP.out.indexed_cram)
        }

        // Graphtyper
        if (params.run_graphtyper) {
            GRAPHTYPER(WEAVER_REMAP.out.sample_name, indexed_genome, WEAVER_REMAP.out.indexed_cram)
        }
    }
}
