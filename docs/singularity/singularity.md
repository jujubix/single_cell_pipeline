# single cell pipeline on a  cluster 

NOTE: the following steps require singularity and access to dockerhub from the compute nodes.


#### input yaml

```
singularity pull docker://docker.io/singlecellpipeline/single_cell_pipeline:v0.4.0
```

### Download the reference data 

The pipeline reference data is available at the following locations:

*Juno cluster at MSKCC:*

```
/juno/work/shah/reference/singlecellpipeline
```

*Shahlab cluster at GSC:*
```
/shahlab/pipelines/reference/singlecellpipeline
```

*Azure:*
```
https://singlecellreference.blob.core.windows.net/refdata
```

NOTE: Access to data on azure is restricted to the shahlab. To request access to the data please contact Diljot Grewal <grewald@mskcc.org.>


### create the context configuration file

The context configuration file contains details about the docker container registry and the directories that need to be mounted inside the singularity containers. Optionally you can also define some default job parameters. The name_match can be used to selectively apply the parameters to the jobs based on their name.

The following example points the pipeline to the dockerhub container registry. The mounts section lists the directories that must be accessible to the pipeline. This includes the directories that contain inputs, output, reference data and input files such as the yaml inputs.
the context section of the yaml snippet below sets the walltime and number of times the job will be retried. The walltime parameter will be set when we launch the pipeline later. 

```
singularity:
    server: 'docker.io'
    username: null
    password: null
    local_cache: '/home/runner/singularity/cache'
    singularity_exe: 'singularity'
    mounts:
      juno: /juno/work/shah/pipelinedir
      reference: /juno/work/shah/reference
      common: /common
context:
  alljobs:
    name_match: '*'
    ctx:
      walltime: '4:00'
      walltime_num_retry: 5
      walltime_retry_increment: '48:00'

```

### launch the pipeline

write the following to the a file:

```
export PATH=/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin:$PATH

single_cell qc --input_yaml /path/to/input.yaml --library_id A97318A --maxjobs 100 \
--sentinel_only  --context_config context.yaml --loglevel DEBUG \
--alignment_output results/alignment --hmmcopy_output results/hmmcopy \
--annotation_output results/annotation --tmpdir temp/temp/QC \
--pipelinedir temp/pipeline/QC  --submit lsf \
--nativespec ' -n {ncpus} -W {walltime} -R "rusage[mem={mem}]span[ptile={ncpus}]select[type==CentOS7]"' \
--config_override '{"refdir": "/path/to/reference/data/dir"}' 
```

launch the pipeline:

```
singularity run --bind /common --bind /juno/work  docker://docker.io/singlecellpipeline/single_cell_pipeline:v0.4.0 sh /path/to/shell/script/from/previous/step
```

The `--bind /common` will mount the `/common` directory inside the singularity. The PATH environment variable must also be set to point to the location of LSF binaries. This will make the commands such as `bsub`, `bjobs` and `bhosts` available to the pipeline. This path will depend on the singularity location. 


This command runs the QC portion of the single cell pipeline. The arguments are:

1. input_yaml: input yaml file containing paths  to fastqs and bam files per cell. 
2. library_id: ID string
3. maxjobs: maximum number of jobs that can be submitted to the scheduler at any given time.
4. sentinel_only: stores file metadata in a sqlite database. pypeliner will look at the timestamps of files when running if this flag is not set which can put a lot of pressure on GPFS storages.
5. context_config: yaml configuration file. see above for details
6. loglevel: DEBUG to print DEBUG level logs
7. alignment_output: takes in a directory path. pipeline will try to align the data if this path is specified
8. hmmcopy_output: takes in a directory path. pipeline will try to run copynumber calling if this path is specified
9. annotation_output: takes in a directory path. pipeline will try to run postprocessing and annotation if this path is specified
10. tmpdir: takes in a directory path, stores intermediate files here
11. pipelinedir: takes in a directory path, stores job level logs and related files here.
12. submit: lsf for LSF clusters, asyncqsub for SGE and local to run it locally on the node
13. nativespec: job specification for the scheduler.  this will depend on the cluster, the settings provided here are for the MSKCC's juno cluster. The values enclosed in brackets will be replaced with the job requirements during job submission. The default values are encoded in the pipeline and can be overridden in the context_config section by the user.
14. config_override: json string that can be used to set pipeline parameters. In this example it is being used to specify the location of reference data that we downloaded in the previous steps.

