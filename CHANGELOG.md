# Change Log

### v0.2.26

##### added
* added fastq screen
	* runs fastqscreen with `--tag`
	* all downstream analysis is run on the tagged data.
	* bam headers contain required information for parsing fastq screen tag
	* by default, pipeline removes all reads that belong to another organism
	* generates a detailed table and adds summary metrics to alignment table 
	* more details at [organism filter](docs/description/organism_filter.md)
* added salmon reference to images.
* added conda package for corrupt tree. updated docker container to use the conda package
* added newick support to heatmap
* added cell order based on corrupt tree to output
* added this changelog
##### changes:
* hmmcopy segments plots have a global max for ylim per run (library)
* standardized page size for corrupt tree output, annotated each page.
* replaced yaml.load with yaml.safe_load
* replaced nan values in QC html with 0
* removed biobloom
* destruct can now handle empty/small fastq files.
* fixed strelka filename issue (missing _)
* refactor alignment workflow
* added a tarball output with all hmmcopy outputs except autoploidy (multipliers 1-6)
* merged all picard based metrics into a single tarball
##### bugs:
* fixed missing header issue with destruct outputs

### v0.2.25
##### added
* Added Corrupt Tree
##### changes
* Reorganized QC pipeline outputs
* updated to newest biobloom container (v0.0.2). biobloom container now runs as root user.
* QC html doesn't require reference GC curve data
* give more memory to biobloom
* load input yaml with safe_load
##### bugs
* bugfix: fixed a merge issue with trim galore running script.

### v0.2.24
##### added
* * Added Cell Cycle Classifier
##### changes
* disable biobloom by default

### v0.2.23
##### bugs
* bugfix: Destruct was not tagging reads with cell ids
##### changes
* removed reference fasta index from github

### v0.2.22
##### changes
* lumpy can now handle empty bams
### v0.2.21
##### added
* added Html QC output
* added biobloom
* added a single 'QC' command to run alignment and hmmcopy
##### changes
* merged the alignment and metrics workflows.
* remove hmmcopy multipliers, only use autoploidy downstream
* removed option to specify multiple hmmcopy parameter sets
##### bugs
* bugfix: issue with automatic dtype detection in csv yaml files.

### v0.2.20
##### changes
* updated input yaml format for pseudowgs. The normal section now follows same schema as tumour (with sample and library id).

### v0.2.19
##### bugs
* bug: missing header in allele_counts file.

### v0.2.18
##### added
* added travis build.
* added smarter dtype merging for csv files.
##### changes
* updated conda recipe
* updated destruct output format from h5 to csv
* fixed destruct to generate counts from filtered output to remove normal reads
* optimized breakpoint calling, normal preprocessing runs only once per run.
* cleaned up raw_dir in output folder
* updated to conda based hmmcopy and mutationseq containers
* updated to latest version of lumpy with correct bed output
* now supports multiple libraries per normal in pseudowgs

### v0.2.17
##### changes
* updated lumpy bed file parsing.
* changed lumpy output file format from h5 to csv.

### v0.2.16
##### changes
* added mutationseq parameters to config. now users can override default settings.

### v0.2.15
##### added
* added parallelization over libraries in pseudowgs
##### bugs
* fixed an issue with read tagging that caused int overflow in bowtie
* some pickling issues due to python compatibility updates in pypeliner
* fixes in csv and yaml generation code

### v0.2.14
##### changes
* refactored lumpy workflow, only run normal preprocessing once per run
* merges in destruct require more disk space
* destruct: read indexes are now unique int
* destuct: reindex both reads in a single job to reduce number of jobs
* destruct: prepocess normal once per run
* faster csv file concatenation
* updated batch config to match pypeliner v0.5.6. now pool selection also accounts for disk usage.
	* Each pool will have available disk space. jobs will be scheduled in a pool based on requirements. production will have smaller disk in standard pool to save on costs.
* classifier now supports csv inputs
##### bugs
* bwa couldnt parse readgroup when not running in docker

### v0.2.13
##### changes
* split and merge bams only when running snv calling
* refactored to make main workflow calling functions standalone subworkflows
* revamped destruct workflow for pseudobulk

### v0.2.12
##### changes
##### changes
* switched to gzipped csv from H5 due to compatibility issues
* order IGV segs file by the clustering order, filter on quality

### v0.2.11
##### added
* added flags to only run parts of pseudowgs workflow
##### bugs
* fixed issue with infer haps where some parameters werent specified correctly.

### v0.2.10
##### changes
* destruct and remixt containers now use the same versioning as single cell pipeline
* lumpy accepts normal cells
* refactor: haplotype calling workflows
* all psudo wgs commands use the same input format as multi sample pseudo bulk.
* bam merge now supports merging larger number of files.

### v0.2.9
##### added
* destruct now supports list of cells as normal
* separate pools based on disk sizes.

##### changes
* separate docker container for destruct
* haplotype calling supports list of cells as normal
* parallel runs support more cells now.
##### bugs
* bug: fixed an issue that caused low mappability mask to disappear in the heatmap

### v0.2.8
##### changes
* replaced python based multiprocessing with gnu parallel
##### bugs
* bug: remixt path fixed in config, mkdir doesnt cause failures anymore in batch vm startup

### v0.2.7
##### added
* added trim galore container
* added a flag to switch disk to 1TB for all batch nodes 
* added a flag to specify whether to trim the fastqs. The flag overrides the sequencer based trimming logic.
* row, column, cell_call and experimental condition can now be null
* switched to gnu parallel for parallel on node runs
* feature: pools are now chosen automatically
##### changes
* updated readgroup string.
* renamed total_mapped_reads column in hmmcopy to total_mapped_reads_hmmcopy to avoid clashed with column of same name in alignment metrics
* snv calling: allow overlaps in vcf files
* remove meta yaml file
* moved autoploidy segment plot to top of page
* added option to launch the pipeline with docker by just adding `--run_with_docker`
* h5 dtype casting uses less memory now
* alignment metrics plot: now faster, plots atmost 1000 cells per page. extra cells overflow onto to the next page.
* support for pypeliner auto detect batch pool
* updated docker 
* updated vcfutils to use pypeliner to handle vcf index files.
* subworkflow resolution now runs in a docker container on compute nodes.
* switched from warnings to logging. the logs from compute now gets reported in main pypeliner log file.
* pipeline now uses OS disk in azure to store temporary files
* updated docker configuration changes in pypeliner. the container doesnt require the prefix anymore.
* added info.yaml file with some metadata per run
* merged andrew's pseudobulk changes.
##### bugs
* fix: strelka uses chromosome size instead of genome size

### v0.2.6
##### changes
* now supports plain text fastq files

### v0.2.5
##### added
* added: heatmap and boxplots for cell quality score
##### changes
* switching to smaller 256GB disks

### v0.2.4
##### changes
* use `table` format in h5 files
##### bugs
* fix: non unique categories error fixed by specifying categories at initialization

### v0.2.3
##### changes
* properly deletes file on batch node after task completes
* cell_id is now a categorical in output
* supports `.fq` and `.fq.gz` fastq file extensions
* heatmap is generated even if all cells are `nan`
* default for non-azure environments is not singularity anymore.

### v0.2.2
##### added
* added docker container info to info yaml files
* added info yaml files 
##### changes
* interprets `?` in picard tools output as 0
* casts all columns in h5 to their correct dtypes.

### v0.2.1
##### added
* added multisample pseudobulk
##### changes
* divided alignment into 2 separate workfloes
* replaced segments and bias pdf with per cell plots with a tarball of png files.

### v0.2.0
##### added
* added singularity support
* added docker support for whole genome
* added docker support for alignment and hmmcopy
* added test data set
##### changes
* project wide refactor
* switched to image with a single disk. removed startup mount commands.
* pseudowgs: added infoer Haplotype code
* clip copynumber in hmmcopy plots to 40

### v0.1.5
##### added
* added LTM
* added haploid poison to hmmcopy
##### changes
* yaml files use block style format
* added cell quality classifier
* bwa-mem is now the default aligner. bwa-aln is also supported
##### bugs
* fix: handles empty segments in hmmcopy segment plot

### v0.1.4
##### added
* pseudowgs: added option to merge bams 
* pseudowgs: added option to split bam by reads (pairs next to each other)
##### changes
* smaller segments plots file size, consistent colormap across plots
* added mask for low mappability regions in heatmap
* one pdf file for segments and bias plots per row of cells 

### v0.1.3
##### added
* added VM image URI and SKU to batch yaml file
* hmmcopy: added autoploidy 
* added titan to pseudowgs
* hmmcopy can now run independently from alignment
##### changes
* rename pick_met to cell_call and condition to experimental_condition
* chooses VM image based on the pipeline version
* cleanly exit hmmcopy script if data is not enough/missing
* aneufinder,alignment, hmmcopy output is now h5
##### bugs
* fixed reouding issue in autoscale formula (there is no round method)

### v0.1.2
##### changes
* streka now runs over split bam files
* can now save split bams using a template specified at run time

### v0.1.1
##### changes
* switch to 1-based state in hmmcopy

### v0.1.0
##### added
* added filtering to copynumber heatmap
* added classifier
* each lane now has a sequencing centre 
##### changes
* metrics heatmaps are not restricted to 72*72
* exposed all hmmcopy params to config file
* modal correction can run on empty datasets
