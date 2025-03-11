## nf-weaver

Nextflow workflow for mapping reads with weaver. The workflow can either map unmapped reads (FASTQ data) or remap previosly mapped reads (BAM/CRAM data) file.

### Usage

Run workflow with

```
nextflow run \
  --map <map.tsv>            # Data to map onto a pangenome with weaver
  --remap <remap.tsv>        # Data to remap onto a pangenome with weaver
  --graph <graph.gfa>        # rGFA or gzipped rGFA
  --genome <genome.fa>       # Indexed FASTA
  --graph_wmi <weaver_index> # Path to the weaver index. If not set, it will be created.
  main.nf
```

or specify these parameters in `conf/params.config`. The `<map.tsv>` should a be a headerless and have two columns: sampleID and fastq path The `<remap.tsv>` should be a headerless and tab separated file with columns: sampleID, sam/bam/cram paths, and list of lanes to remap ("." if all lanes should be remapped). Each line should have a unique sampleID and the sam paths are separeted by whitespace (' ') and are expected to a be in SAM.gz/BAM/CRAM format. When selected lanes are remapped only, they should be also separated by whitespace (' '). The `<graph.gfa>` is the reference pangenome graph in rGFA format and may be gzipped. The `<genome.fa>` is an indexed FASTA file.


#### Reference cache parameters

For CRAM encoding/decoding tasks it is recommended to specify the location of a reference cache. These are specified using these optional parameters:

```
  --ref_cache_tar_gz <tar.gz> # Provide archive containing the reference cache for CRAM encoding/decoding
  --ref_cache <tar.gz>        # Provide path to the reference cache for CRAM encoding/decoding
```

Either the --ref_cache_tar_gz parameter, --ref_cache parameter or `REF_CACHE` environment variable should be set for CRAM encoding/decoding to work properly.


#### Weaver index creation parameters

If the workflow is run without providing a weaver index input, it will be created first. These optional parameters will affect how the index is created:

```
  --weaver_k <k> # Use this `k`-mer size when creating a weaver index.
  --weaver_w <w> # Find minimizer among `w` many k-mers when creating a weaver index.
  --vcf <vcf.gz> # Use variants from this VCF.gz when creating a weaver index.
```

#### Debugging

If case the workflow consistently fails, you might want to set `params.scratch = false` (in your custom config or `conf/params.config`) to access all logs.

## License

MIT
