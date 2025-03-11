/*!
 * weaver index module
 *
 * Indexes a reference pangenome graph (rGFA) with input small variants in VCF.gz.
 */

params.weaver_index_label = "weaver_index"
params.weaver_index_publish_dir = "${launchDir}/${params.weaver_index_label}"
params.weaver_index_publish_mode = "symlink"

//! Process for running 'weaver index'
process WEAVER_INDEX {
    publishDir "${params.weaver_index_publish_dir}", mode: "${params.weaver_index_publish_mode}"
    label "${params.weaver_index_label}"

    input:
    //! Reference pangenome graph (rGFA)
    val graph

    //! k-mer and minimizer size
    val k

    //! Number of k-mers in a window
    val w

    output:
    tuple path("weaver_indexed.gfa"), path("weaver_indexed.gfa.wmi"), emit: indexed_graph

    script:
    """
    weaver=${params.weaver_path}

    cp -a ${graph} weaver_indexed.gfa

    if [[ -e ${graph}.alts ]]; then
      ln -s ${graph}.alts weaver_indexed.gfa.alts
    fi

    if [[ ! -e "${params.graph_wmi}" ]]; then
      \${weaver} index weaver_indexed.gfa \
        --vcf=${params.vcf} \
        --k=${params.weaver_k} \
        --w=${params.weaver_w} \
        --threads=${task.cpus} \
        --log=weaver_index.log \
        --verbose
    else
      echo "INFO: weaver index skipped"
      ln -s ${params.graph_wmi} weaver_indexed.gfa.wmi
    fi
    """
}
