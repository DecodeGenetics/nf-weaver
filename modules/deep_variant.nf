/**
 * DeepVariant genotyping process for weaver mapped data.
 */

params.deep_variant_label = "deep_variant"
params.deep_variant_publish_dir = "${launchDir}/${params.deep_variant_label}"
params.deep_variant_publish_mode = "symlink"

process DEEP_VARIANT {
    label "${params.deep_variant_label}"
    publishDir("${params.deep_variant_publish_dir}", mode: "${params.deep_variant_publish_mode}")
    tag "${sample_name}"

    input:
    val(sample_name)
    tuple path(genome), path(genome_fai)
    tuple path(cram), path(cram_crai)

    output:
    //! DeepVariant genotype calls
    tuple path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit: indexed_vcf

    //! Add extra deep variant options
    def extra_deep_variant_opts = (params.customized_dv_model != "") ? "--customized_model=${params.customized_dv_model}" : ""

    script:
    """
    bash /opt/deepvariant/bin/run_deepvariant ${extra_deep_variant_opts} \
      --model_type=WGS \
      --ref=${genome} \
      --reads=${cram} \
      --use_keras_model=true \
      --output_vcf=${sample_name}.vcf.gz \
      --num_shards=${task.cpus} \
      --sample_name=${sample_name} \
      --make_examples_extra_args="min_mapping_quality=0,normalize_reads=true"
    """
}
