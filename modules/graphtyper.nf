/**
 * graphtyper genotyping process for weaver mapped data.
 */

params.graphtyper_label = "graphtyper"
params.graphtyper_publish_dir = "${launchDir}/${params.graphtyper_label}"
params.graphtyper_publish_mode = "symlink"

process GRAPHTYPER {
    label "${params.graphtyper_label}"
    publishDir "${params.graphtyper_publish_dir}"
    tag "${sample_name}"

    input:
    val(sample_name)
    tuple path(genome), path(genome_fai)
    tuple path(cram), path(cram_crai)

    output:
    //! graphtyper genotype calls in VCF.gz and its index
    tuple path("${sample_name}.vcf.gz"), path("${sample_name}.vcf.gz.tbi"), emit: indexed_vcf

    def extra_deep_variant_opts = (params.customized_dv_model != "") ? "--customized_model=${params.customized_dv_model}" : ""

    script:
    if (params.bed_file != "")
        error "BED file mode is not yet implemented. Talk to a developer."
    else
        """
        for c in {1..22} X
        do
          b=1
          e=\$(grep -m1 -w "chr\${c}" ${genome_fai} | cut -f2)

          while [[ \$b -lt \$e ]]; do
            be=\$((b + 4999999 > e ? e : b + 4999999))
            echo "graphtyper genotype ${genome} --sam=${cram} --threads=1 --region chr\${c}:\${b}-\${be}"
            b=\$((be + 1))
          done
        done | tr '\\n' '\\0' | xargs -0 -n1 -P${task.cpus} sh -c

        for c in {1..22} X; do
          if [[ -e "results/chr\${c}" ]]; then
            find "results/chr\${c}" -name "*.vcf.gz" | sort
          fi
        done > ./vcf_file_list

        if true; then
          bcftools view --no-version --header-only "\$(head -n1 vcf_file_list)"

          while read -r vcf ; do
            bgzip -dc "\${vcf}" | grep -v ^#
          done < ./vcf_file_list
        fi | bgzip -@${task.cpus} -c > ${sample_name}.vcf.gz

        tabix -p vcf ${sample_name}.vcf.gz
        """
}
