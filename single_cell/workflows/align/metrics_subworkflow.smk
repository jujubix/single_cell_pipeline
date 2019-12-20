import os
from single_cell.workflows.align.dtypes import dtypes

rule get_duplication_wgs_flagstat_metrics:
    input:
        input_bam = os.path.join(bams_dir, '{cell_id}.bam')
    params:
        ref_genome = ref_genome,
        picard_wgs_params = alignment_config['picard_wgs_params'],
        picard_docker = alignment_config.get('docker', {}).get('picard', None)
    output:
        markdups_bam = temp('{cell_id}_temp_markdup_bam.bam'),
        markdups_metrics = temp('{cell_id}_markdups_metrics.txt'),
        tempdir = temp(directory('{cell_id}_tempdir_markdups')),
        wgs_metrics = temp('{cell_id}_wgs_metrics.txt')
    run:
        single_cell.workflows.align.tasks.picard_wgs_dup(
        input.input_bam, output.markdups_bam, output.markdups_metrics, output.tempdir,
        params.ref_genome, output.wgs_metrics, params.picard_wgs_params,
        params.picard_docker
        )

rule bam_collect_gc_insert_metrics:
    input:
        input_bam = os.path.join(bams_dir, '{cell_id}.bam'),
    params:
        ref_genome = ref_genome,
        picard_docker = alignment_config.get('docker', {}).get('picard', None),
        samtools_docker = alignment_config.get('docker', {}).get('samtools', None)
    output:
        gc_metrics = temp('{cell_id}_gc_metrics.txt'),
        gc_metrics_summary = temp('{cell_id}_gc_metrics_summary.txt'),
        gc_metrics_pdf = temp('{cell_id}_gc_metrics.pdf'),
        tempdir = temp(directory('{cell_id}_gc_tempdir')),
        flagstat_metrics = temp('{cell_id}_flagstat_metrics.txt'),
        insert_metrics = temp('{cell_id}_insert_metrics.txt'),
        insert_pdf = temp('{cell_id}_insert_metrics.pdf')
    resources:
        mem = alignment_config.get('memory', {}).get('med', 6),
        ncpus = 1
    run:
        single_cell.workflows.align.tasks.picard_insert_gc_flagstat(
        input.input_bam, params.ref_genome, output.gc_metrics,
        output.gc_metrics_summary, output.gc_metrics_pdf,
        output.tempdir, output.flagstat_metrics, output.insert_metrics,
        output.insert_pdf, params.picard_docker, params.samtools_docker
        )

rule collect_gc_metrics:
    input:
        expand('{cell_id}_gc_metrics.txt', cell_id = cell_id)
    output:
        outfile = alignment_files['gc_metrics_csv'],
        outfile_yaml = '{}.yaml'.format(alignment_files['gc_metrics_csv']),
        tempdir = temp(directory("temp_gc"))
    resources:
        mem = alignment_config.get('memory', {}).get('med', 6),
        ncpus = 1
    run:
        input_dict = {cid: f'{cid}_gc_metrics.txt' for cid in cell_id}
        single_cell.workflows.align.tasks.collect_gc(input_dict, output.outfile, output.tempdir)

rule collect_metrics:
    input:
        flagstat_metrics = expand('{cell_id}_flagstat_metrics.txt', cell_id = cell_id),
        markdups_metrics = expand('{cell_id}_markdups_metrics.txt', cell_id = cell_id),
        insert_metrics = expand('{cell_id}_insert_metrics.txt', cell_id = cell_id),
        wgs_metrics = expand('{cell_id}_wgs_metrics.txt', cell_id = cell_id)
    output:
        tempdir = temp(directory("tempdir_collect_metrics")),
        merged_metrics = temp("alignment_metrics.csv.gz"),
        merged_metrics_yaml = temp("alignment_metrics.csv.gz.yaml")
    resources:
        mem = alignment_config.get('memory', {}).get('med', 6),
        ncpus = 1
    run:
        flagstat_metrics = {cid: f'{cid}_flagstat_metrics.txt' for cid in cell_id}
        markdups_metrics = {cid: f'{cid}_markdups_metrics.txt' for cid in cell_id}
        insert_metrics = {cid: f'{cid}_insert_metrics.txt' for cid in cell_id}
        wgs_metrics = {cid: f'{cid}_wgs_metrics.txt' for cid in cell_id}
        single_cell.workflows.align.tasks.collect_metrics(
            flagstat_metrics, markdups_metrics, insert_metrics,
            wgs_metrics, output.tempdir, output.merged_metrics
            )

rule annotate_metrics:
    input:
        infile = "alignment_metrics.csv.gz"
    params:
        annotation_data = sampleinfo,
        dtypes = dtypes()['metrics']
    output:
        outfile = temp('alignment_metrics_annotated.csv.gz'),
        outfile_yaml = temp('alignment_metrics_annotated.csv.gz.yaml')
    resources:
        mem = alignment_config.get('memory', {}).get('med', 6),
        ncpus = 1
    run:
        single_cell.utils.csvutils.annotate_csv(
            input.infile, params.annotation_data, output.outfile, dtypes=params.dtypes
            )

rule add_fastqscreen_metrics:
    input:
        in_filenames = [
            'alignment_metrics_annotated.csv.gz',
            'organism_summary_count_per_cell.csv'
        ]
    output:
        out_filename = temp('alignment_metrics.csv'),
        out_filename_yaml = temp('alignment_metrics.csv.yaml')
    params:
        how = 'outer',
        on = ['cell_id'],
        dtypes = dtypes()['metrics']
    resources:
        mem = alignment_config.get('memory', {}).get('med', 6),
        ncpus = 1
    run:
        single_cell.utils.csvutils.merge_csv(
            input.in_filenames, output.out_filename,
            params.how, params.on, dtypes=params.dtypes
            )

