/**
 * weaver map module. Used for mapping reads from input FASTQ file(s) using weaver.
 *
 * By default the workflow will run in "lossy" mode, which will remove read duplicates, remove unplaced reads,
 * and use the tool 'crumble' on the output CRAM file. Use the input parameter '--lossless' if you want the mapping
 * to be lossless.
 */

params.weaver_map_label = "weaver_map"
params.weaver_map_publish_dir = "${launchDir}/${params.weaver_map_label}"
params.weaver_map_publish_mode = "symlink"

//! Process for remapping FASTQ file(s) using weaver.
process WEAVER_MAP {
    publishDir("${params.weaver_map_publish_dir}", mode: "${params.weaver_map_publish_mode}")
    label "${params.weaver_map_label}"
    tag "${sample_name}"

    input:
    //! Reference pangenome and its index. Typically in GFA or GFA.gz, but can also be FASTA or FASTA.gz
    tuple path(graph), path(graph_wmi)

    //! FASTA genome. Required to be indexed. Usually this is the FASTA that is transformed from an rGFA.
    val genome

    /*!
     * Sample and and its input path(s)
     *   sample_name: Name of sample, is used to name the output file.
     *   '\t'       : Tab character in between.
     *   fastq(s)   : One or two FASTQ files to map, separated by a whitespace (' ') when two.
     */
    tuple val(sample_name), val(fastqs)

    //! Remove duplicates in lossy mode
    def extra_samtools_markdup_opts = params.lossless ? "" : "-r"

    //! Bin quality values in lossy mode
    def extra_weaver_map_opts = params.lossless ? "" : "--rta3-quals"

    //! Option to use when crumbling a CRAM file
    def crumble_cram_opts = params.crumble_cram_opts != "" ?
        params.crumble_cram_opts :
        "-1 -OCRAM,version=3.0,lossy_names=1,seqs_per_slice=25000"

    //! Set a specific read group header line. If empty, an automatically one will be created.
    def read_group_header_line = params.read_group_header_line ? params.read_group_header_line : ""

    //! Set to enable marking duplicates (default true).
    def is_marking_duplicates = params.is_marking_duplicates ? "true" : "false"

    output:
    //! Sample name identifier.
    val sample_name, emit: sample_name

    //! CRAM file mapped with weaver and its index
    tuple path("${sample_name}.cram"), path("${sample_name}.cram.crai"), emit: indexed_cram

    //! Weaver log output
    path "${sample_name}.weaver_map.log", emit: weaver_log

    script:
    """
    weaver=${params.weaver_path} # For development, it is possible to set a custom weaver path

    pwd=\$(pwd)
    echo "Hostname:  \$(hostname -s)"
    echo "Timestamp: \$(date)"
    echo "Dir:       \${pwd}"
    echo "REF_CACHE: \${REF_CACHE}"
    echo "weaver:    \${weaver}"

    ## Create temporary directories
    mkdir markdup # Temporary directory for samtools markdup
    mkdir crumble # Temporary directory for crumble

    ## FASTQ arguments
    fastq1=\$(echo ${fastqs} | awk -v IFS=' ' 'END{print \$1}')
    fastq2=\$(echo ${fastqs} | awk -v IFS=' ' 'END{print \$2}')

    if [[ -n \${fastq2} ]]
    then
      weaver_map_fastq_opts=("\${fastq1}" "--fq2=\${fastq2}")
    else
      weaver_map_fastq_opts=("\${fastq1}")
    fi

    ## Read group options
    if [[ -n "${read_group_header_line}" ]]
    then
      read_group_header_line=\'${read_group_header_line}\'
    else
      read_group_header_line='@RG\\tID:1\\tSM:${sample_name}'
    fi

    if ${is_marking_duplicates}; then
      ## weaver map > mark duplicates > write CRAM
      \${weaver} map ${graph} \${weaver_map_fastq_opts} --threads=${task.cpus} --log=${sample_name}.weaver_map.log ${extra_weaver_map_opts} --read-group-header-line=\${read_group_header_line} --verbose |
        samtools markdup -l200 -Tmarkdup/ -u ${extra_samtools_markdup_opts} - - |
        samtools view --remove-tag ms,MC -OCRAM,version=3.0,seqs_per_slice=25000,reference=${genome} -o weaver.cram -
    else
      ## weaver map > write CRAM
      \${weaver} map ${graph} \${weaver_map_fastq_opts} --threads=${task.cpus} --log=${sample_name}.weaver_map.log ${extra_weaver_map_opts} --read-group-header-line=\${read_group_header_line} --verbose |
        samtools view --remove-tag ms,MC -OCRAM,version=3.0,seqs_per_slice=25000,reference=${genome} -o weaver.cram -
    fi

    if ${params.lossless}
    then
      ## Skip crumble for lossless output
      mv weaver.cram ${sample_name}.cram
    else
      ## Use Crumble if the output can be lossy
      ## Index the output CRAM file
      samtools index weaver.cram

      ## Find all contigs that have some sequences mapped
      gzip -dc weaver.cram.crai | cut -f1 | uniq > contigs_with_seqs.txt

      ## Run Crumble in parallel for each input contig with any reads (largest contigs first)
      samtools view --no-PG -H weaver.cram |
        awk '\$1 == "@SQ"{
        sn=""; ln=0;
        for (i = 2; i <= NF; i++){
          if      (substr(\$i,1,3) == "SN:") { sn = substr(\$i,4) }
          else if (substr(\$i,1,3) == "LN:") { ln = substr(\$i,4) }
        }
        print sn" "ln" "NR-1
      }' |
        IFS=" " sort -k2,2nr -k1,1 |
        awk -v FS=' ' 'NR==FNR{m[\$1]=1} NR!=FNR && (\$3 in m){printf "%07d %s\\n", \$3, \$1}' contigs_with_seqs.txt - |
        xargs -P${task.cpus} -n2 bash -c 'crumble ${crumble_cram_opts},reference=${genome} -r "\${2}" weaver.cram crumble/\${1}.cram' bash

      ## Combine all region and write final output
      find crumble/ -name "*.cram" | LC_ALL=C sort | xargs samtools cat -o ${sample_name}.cram
    fi

    ## Clean temporary directories in the background
    rm -rf -- \$tmp &

    ## Index the output CRAM
    samtools index ${sample_name}.cram

    ## Wait my love
    wait
    """
}
