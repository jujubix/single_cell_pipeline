'''
Created on Jun 6, 2018

@author: dgrewal
'''
import yaml
import config_reference
import collections

def override_config(config, override):
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        return d

    if not override:
        return config

    cfg = update(config, override)

    return cfg


def get_config_params(override=None):
    input_params = {
        "cluster": "azure", "aligner": "bwa-mem",
        "reference": "grch37", "smoothing_function": "modal",
        "bin_size": 500000, "copynumber_bin_size": 1000,
        'memory': {'high': 18, 'med': 6, 'low': 2}
    }

    input_params = override_config(input_params, override)

    return input_params

def write_config(params, filepath):
    with open(filepath, 'w') as outputfile:
        yaml.safe_dump(params, outputfile, default_flow_style=False)


def get_hmmcopy_params(cluster, reference, binsize, smoothing_function):
    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']
    docker_containers = {
        'single_cell_pipeline': docker_containers['single_cell_pipeline'],
        'hmmcopy': docker_containers['hmmcopy']
    }

    params = {
        'multipliers': [1, 2, 3, 4, 5, 6],
        'map_cutoff': 0.9,
        'bin_size': binsize,
        'e': 0.999999,
        'eta': 50000,
        'g': 3,
        'lambda': 20,
        'min_mqual': 20,
        'nu': 2.1,
        'num_states': 12,
        's': 1,
        'strength': 1000,
        'kappa': '100,100,700,100,25,25,25,25,25,25,25,25',
        'm': '0,1,2,3,4,5,6,7,8,9,10,11',
        'mu': '0,1,2,3,4,5,6,7,8,9,10,11',
        'smoothing_function': smoothing_function,
        'exclude_list': referencedata['exclude_list'],
        'gc_wig_file': referencedata['gc_wig_file'][binsize],
        'map_wig_file': referencedata['map_wig_file'][binsize],
        'classifier_training_data': referencedata['classifier_training_data'],
        'chromosomes': referencedata['chromosomes'],
        'ref_genome': referencedata['ref_genome'],
        'docker': docker_containers,
        'memory': {'med': 6},
        'good_cells': [['median_hmmcopy_reads_per_bin', 'ge', 50]]
    }

    return {"hmmcopy": {"autoploidy": params}}


def get_align_params(cluster, reference, binsize, smoothing_function, aligner):
    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']
    docker_containers = {
        'single_cell_pipeline': docker_containers['single_cell_pipeline'],
        'fastqc': docker_containers['fastqc'],
        'samtools': docker_containers['samtools'],
        'bwa': docker_containers['bwa'],
        'picard': docker_containers['picard'],
    }

    params = {
        'ref_genome': referencedata['ref_genome'],
        'docker': docker_containers,
        'memory': {'med': 6},
        'aligner': aligner,
        'adapter': 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        'adapter2': 'CTGTCTCTTATACACATCTGACGCTGCCGACGA',
        'picard_wgs_params' : {
            "min_bqual": 20,
            "min_mqual": 20,
            "count_unpaired": False,
        },
        'chromosomes': referencedata['chromosomes'],
        'gc_windows': referencedata['gc_windows']
    }

    return {"alignment": params}


def get_aneufinder_params(cluster, reference):
    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']
    params = {
        'memory': {'med': 6},
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'aneufinder': docker_containers['aneufinder'],
        },
        'chromosomes': referencedata['chromosomes'],
        'ref_genome': referencedata['ref_genome']
    }

    return {'aneufinder': params}


def get_merge_bams_params(cluster, reference):

    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']
    params = {
        'memory': {'med': 6, 'high': 18},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools']
        },
        'ref_genome': referencedata['ref_genome'],
        'split_size': 10000000,
        'chromosomes': referencedata['chromosomes'],
    }
    return {'merge_bams': params}


def get_split_bam_params(cluster, reference):

    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']

    params = {
        'memory': {'low':4, 'med': 6, 'high': 18},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools']
        },
        'ref_genome': referencedata['ref_genome'],
        'split_size': 10000000,
        'chromosomes': referencedata['chromosomes'],
        'one_split_job': False
    }

    return {'split_bam': params}


def get_germline_calling_params(cluster, reference):

    if cluster == "azure":
        referencedata = config_reference.reference_data_azure(reference)
    else:
        referencedata = config_reference.reference_data_shahlab(reference)

    docker_containers = config_reference.containers()['docker']

    params = {
        'memory': {'low': 4, 'med': 6, 'high': 18},
        'max_cores': 8,
        'docker': {
            'single_cell_pipeline': docker_containers['single_cell_pipeline'],
            'samtools': docker_containers['samtools'],
            'vcftools': docker_containers['vcftools'],
            'snpeff': docker_containers['snpeff'],
        },
        'ref_genome': referencedata['ref_genome'],
        'chromosomes': referencedata['chromosomes'],
        'split_size': 10000000,
        'databases':{
            'mappability':{
               'url': 'http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/release3'
                      '/wgEncodeCrgMapabilityAlign50mer.bigWig',
               'local_path': referencedata['databases']['mappability']['local_path'],
            },
            'snpeff': {"db": 'GRCh37.75'},
        },
    }

    return {'germline_calling': params}



def get_singlecell_pipeline_config(config_params, override=None):
    reference = config_params["reference"]
    cluster = config_params["cluster"]

    params = {}

    params.update(
        get_hmmcopy_params(
            cluster, reference, config_params["bin_size"],
            config_params["smoothing_function"]
        )
    )

    params.update(
        get_align_params(
            cluster, reference, config_params["bin_size"],
            config_params["smoothing_function"],
            config_params['aligner'],

        )
    )

    params.update(get_aneufinder_params(cluster, reference))

    params.update(get_merge_bams_params(cluster, reference))

    params.update(get_split_bam_params(cluster, reference))

    params.update(get_germline_calling_params(cluster, reference))

    params = override_config(params, override)

    return params


