'''
Created on Feb 22, 2018

@author: dgrewal
'''
import os
import pypeliner
import pypeliner.managed as mgd
from workflows import pseudo_wgs 
from single_cell.utils import helpers

def pseudo_wgs_workflow(workflow, args):

    config = helpers.load_config(args)
    bam_files, _  = helpers.get_bams(args['input_yaml'])
    sampleids = helpers.get_samples(args['input_yaml'])


    wgs_bam_dir = args["merged_wgs"]

    wgs_bam_template = os.path.join(wgs_bam_dir, "{regions}_merged.bam")
    wgs_bai_template = os.path.join(wgs_bam_dir, "{regions}_merged.bam.bai")

    


    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sampleids,
    )

    workflow.transform(
        name="get_regions",
        ctx={'mem': 2, 'num_retry': 3, 'mem_retry_increment': 2, 'pool_id': config['pools']['standard'], 'ncpus':1 },
        func=helpers.get_bam_regions,
        ret=pypeliner.managed.TempOutputObj('regions'),
        args=(
              config["ref_genome"],
              int(10**7),
              config["chromosomes"],
        )
    )

    workflow.subworkflow(
        name='wgs_workflow',
        func=pseudo_wgs.create_wgs_workflow,
        args=(
            mgd.InputFile('bam_markdups', 'sample_id', fnames=bam_files),
            mgd.OutputFile("merged_bam", "regions", axes_origin=[], template=wgs_bam_template),
            mgd.OutputFile("merged_bai", "regions", axes_origin=[], template=wgs_bai_template),
            sampleids,
            config,
            mgd.TempInputObj("regions"),
        )
    )


    return workflow