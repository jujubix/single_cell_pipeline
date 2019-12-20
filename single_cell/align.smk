import os
import re
import sys
from single_cell.utils import inpututils

def get_output_files(outdir, lib):
    data = {
        'alignment_metrics_csv': os.path.join(outdir, '{}_alignment_metrics.csv.gz'.format(lib)),
        'gc_metrics_csv': os.path.join(outdir, '{}_gc_metrics.csv.gz'.format(lib)),
        'fastqc_metrics_csv': os.path.join(outdir, '{}_detailed_fastqscreen_metrics.csv.gz'.format(lib)),
        'plot_metrics_output': os.path.join(outdir, '{}_alignment_metrics.pdf'.format(lib)),
        'alignment_metrics_tar': os.path.join(outdir, '{}_alignment_metrics.tar.gz'.format(lib)),
    }

    return data

alignment_config = inpututils.load_config(config)
alignment_config = alignment_config['alignment']

lib = config["library_id"]
alignment_dir = config["out_dir"]
bams_dir = config["bams_dir"]

sampleinfo = inpututils.get_sample_info(config['input_yaml'])
laneinfo = inpututils.get_lane_info(config['input_yaml'])

cellids = inpututils.get_samples(config['input_yaml'])
fastq1_files, fastq2_files = inpututils.get_fastqs(config['input_yaml'])

alignment_files = get_output_files(alignment_dir, lib)
alignment_meta = os.path.join(alignment_dir, 'metadata.yaml')

bam_files_template = os.path.join(bams_dir, '{cell_id}.bam')
bams_meta = os.path.join(bams_dir, 'metadata.yaml')

lanes = sorted(set([v[1] for v in fastq1_files.keys()]))
cell_id = sorted(set([v[0] for v in fastq1_files.keys()]))

input_yaml_blob = os.path.join(alignment_dir, 'input.yaml')

baseimage = alignment_config.get('docker', {}).get('single_cell_pipeline', None)

chromosomes = alignment_config["chromosomes"]
ref_genome = alignment_config['ref_genome']

include: "workflows/align/Snakefile"

rule all:
    input:
        bam_markdups = expand(os.path.join(bams_dir, '{cell_id}.bam'), cell_id = cell_id),
        alignment_metrics_csv = alignment_files['alignment_metrics_csv'],
        gc_metrics_csv = alignment_files['gc_metrics_csv'],
        fastqc_metrics_csv = alignment_files['fastqc_metrics_csv'],
        plot_metrics_output = alignment_files['plot_metrics_output'],
        alignment_meta = alignment_meta,
        bams_meta = bams_meta

rule generate_meta_files_results:
    params:
        command = sys.argv[0:],
        root_dir = alignment_dir,
        input_yaml_data = inpututils.load_yaml(config['input_yaml']),
        input_yaml = input_yaml_blob,
        metadata = {
                'library_id': lib,
                'cell_ids': cell_id,
                'lane_ids': lanes,
                'type': 'alignment'
            }
    input:
        filepaths=list(alignment_files.values())
    output:
        alignment_meta
    run:
        single_cell.utils.helpers.generate_and_upload_metadata(
            params.command, params.root_dir, input.filepaths, output,
            input_yaml_data = params.input_yaml_data, input_yaml = params.input_yaml,
            metadata = params.metadata
            )

rule generate_meta_files_bams:
    params:
        command = sys.argv[0:],
        root_dir = bams_dir,
        template = (cell_id, os.path.join(bams_dir, '{cell_id}.bam'), 'cell_id'),
        metadata = {
                'library_id': lib,
                'cell_ids': cell_id,
                'lane_ids': lanes,
                'type': 'cellbams'
            }
    input:
        filepaths = expand(os.path.join(bams_dir, '{cell_id}.bam'), cell_id = cell_id)
    output:
        bams_meta
    run:
        single_cell.utils.helpers.generate_and_upload_metadata(
            params.command, params.root_dir, input.filepaths, output,
            template = params.template, metadata = params.metadata
            )
