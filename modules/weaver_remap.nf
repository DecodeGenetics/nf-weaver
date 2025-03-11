/**
 * weaver remap module. For remapping reads from input SAM/BAM/CRAM files using weaver.
 *
 * By default the workflow will run in "lossy" mode, which will remove read duplicates, remove unplaced reads,
 * and use the tool 'crumble' on the output CRAM file. Use the input parameter '--lossless' if you want the mapping
 * to be lossless.
 */

params.weaver_remap_label = "weaver_remap"
params.weaver_remap_publish_dir = "${launchDir}/${params.weaver_remap_label}"
params.weaver_remap_publish_mode = "symlink"

//! Process for remapping a SAM/BAM/CRAM file using weaver.
process WEAVER_REMAP {
    publishDir("${params.weaver_remap_publish_dir}", mode: "${params.weaver_remap_publish_mode}")
    label "${params.weaver_remap_label}"
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
     *   sams       : SAM/BAM/CRAM files to remap, separated by a whitespace (' ').
     *               'sam' may be position sorted. The output will be merged when there are more than one file listed.
     *   lanes      : Lanes to extract for remapping, separated by a whitespace (' '). Set as "." for all lanes.
     */
    tuple val(sample_name), val(sams), val(lanes)

    //! Remove duplicates in lossy mode
    def extra_samtools_markdup_opts = params.lossless ? "" : "-r"

    //! Bin quality values in lossy mode
    def extra_weaver_map_opts = params.lossless ? "" : "--rta3-quals"

    //! Option to use when crumbling a CRAM file
    def crumble_cram_opts = params.crumble_cram_opts != "" ?
        params.crumble_cram_opts :
        "-1 -OCRAM,version=3.0,lossy_names=1,seqs_per_slice=25000"

    def is_overwriting_sample_name_in_rg_string = params.is_overwriting_sample_name_in_rg_string ? "true" : "false"

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
    echo "sams:      ${sams}"
    echo "lanes:     ${lanes}"

    ## Create temporary directories
    mkdir collate # Temporary directory with collated reads
    mkdir markdup # Temporary directory for samtools markdup
    mkdir crumble # Temporary directory for crumble

    samtools_view_lanes_opts=("--no-PG")

    ## Subset lanes
    if [[ "${lanes}" != "." ]]; then
      echo "${lanes}" | tr ' ' '\\n' > lanes.tsv # File containing the lanes to subset on
      samtools_view_lanes_opts+=("-R" "lanes.tsv")
      echo "INFO: Number of lanes to subset to is \$(cat lanes.tsv | wc -l)"
    fi

    ## Get read group header lines
    for sam in ${sams}; do
      samtools view "\${samtools_view_lanes_opts[@]}" -H "\${sam}" | grep '^@RG'
    done > extra_header_lines.tsv

    if ${is_overwriting_sample_name_in_rg_string}; then
      sed -i 's;\tSM:[[:print:]]*;\tSM:${sample_name};' extra_header_lines.tsv
    fi

    echo "INFO: Number of read group header lines are \$(cat extra_header_lines.tsv | wc -l)"

    ## Collate > Convert to FASTQ > weaver map > mark duplicates > write CRAM
    for sam in ${sams}
    do
      ## Collate/group read pairs together
      echo "INFO: Collating \${sam} ..."
      samtools view -u "\${samtools_view_lanes_opts[@]}" "\${sam}" |
        samtools collate -O -u -n128 -Tcollate/ - |
        samtools fastq -0 /dev/null -s /dev/null -T RG - # Print the FASTQ data to stdout
    done |
      \${weaver} map ${graph} - --threads=${task.cpus} --log=${sample_name}.weaver_map.log ${extra_weaver_map_opts} --verbose --extra-header-lines=extra_header_lines.tsv |
      samtools markdup -l200 -Tmarkdup/ --use-read-groups -u ${extra_samtools_markdup_opts} - - |
      samtools view --remove-tag ms,MC -OCRAM,version=3.0,seqs_per_slice=25000,reference=${genome} -o weaver.cram -

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
