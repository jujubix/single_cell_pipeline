"""
Extract metrics table for multiple samples.

TODO:
- Option to also output df_hist for wgs metrics and insert size metrics
- Replace out_file with out_basename, output three files?
"""

from __future__ import division
from collections import OrderedDict

import argparse
import os
import pandas as pd

#=======================================================================================================================
# Read Command Line Input
#=======================================================================================================================
parser = argparse.ArgumentParser()

parser.add_argument('metrics_dir',
                    help='''Path to metrics directory generated by alignment pipeline.''')

parser.add_argument('out_file',
                    help='''Path to .csv file where table output will be written.''')

parser.add_argument('--library_id',
                    help='''Optional identifier string for the library.''')

args = parser.parse_args()

#=======================================================================================================================
# Functions
#=======================================================================================================================

def extract_wgs_metrics(dir):
    ''' Extract coverage depth and breadth '''
    
    metrics = pd.DataFrame()
    
    for file in os.listdir(dir):
        if file.endswith('.rmdups.txt'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep='\t', header=6, skipfooter=1, engine='python')
            
            df_metrics = df.iloc[[0]]
            df_metrics.columns = [x.lower() for x in df_metrics.columns]
            
            df_hist = df.iloc[4:,0:2]
            df_hist.columns = ['coverage', 'count']
            df_hist = df_hist.reset_index(drop=True)
            
            genome_territory = int(df_metrics.ix[0, 'genome_territory'])
            
            coverage_depth = df_metrics.ix[0, 'mean_coverage']
            coverage_breadth = (genome_territory - int(df_hist.ix[0, 'count'])) / genome_territory
            
            sample_metrics = pd.DataFrame(OrderedDict([
                                                       ('sample_id', [sample_id]), 
                                                       ('coverage_depth', [coverage_depth]), 
                                                       ('coverage_breadth', [coverage_breadth])
                                                       ]))
            
            metrics = metrics.append(sample_metrics)
    
    metrics = metrics.reset_index(drop=True)
    
    return metrics

def extract_insert_metrics(dir):
    ''' Extract median and mean insert size '''
    
    metrics = pd.DataFrame()
    
    for file in os.listdir(dir):
        if file.endswith('.rmdups.txt'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep='\t', header=6, skipfooter=1, engine='python')
            
            df_metrics = df.iloc[[0]]
            df_metrics.columns = [x.lower() for x in df_metrics.columns]
            
            df_hist = df.iloc[4:,0:2]
            df_hist.columns = ['size', 'count']
            df_hist = df_hist.reset_index(drop=True)
            
            sample_metrics = pd.DataFrame(OrderedDict([
                                                       ('sample_id', [sample_id]), 
                                                       ('median_insert_size', df_metrics.ix[0, 'median_insert_size']), 
                                                       ('mean_insert_size', df_metrics.ix[0, 'mean_insert_size']),
                                                       ('standard_deviation_insert_size', df_metrics.ix[0, 'standard_deviation'])
                                                       ]))
            
            metrics = metrics.append(sample_metrics)
    
    metrics = metrics.reset_index(drop=True)
        
    return metrics

def extract_duplication_metrics(dir):
    ''' Extract duplication metrics and library size '''
    
    metrics = pd.DataFrame()
    
    for file in os.listdir(dir):
        if file.endswith('.markdups.txt'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep='\t', header=6, skipfooter=1, engine='python')
            
            df.columns = [x.lower() for x in df.columns]
            
            if any(df['library'] == '## HISTOGRAM'):
                hist_index = df['library'][df['library'] == '## HISTOGRAM'].index.tolist()[0]
                
                df_metrics = df[0:hist_index-1]
            else: 
                df_metrics = df[0:-1]
            
            if len(df) > 1:
                df_metrics = df_metrics.drop(['library', 'percent_duplication'], axis=1).sum(axis=0)
                
                df_metrics['percent_duplication'] = (df_metrics['unpaired_read_duplicates'] + 
                                                   ((df_metrics['read_pair_duplicates'] + df_metrics['read_pair_optical_duplicates']) * 2)) / \
                                                    (df_metrics['unpaired_reads_examined'] + (df_metrics['read_pairs_examined'] * 2))
                
                df_metrics = pd.DataFrame(df_metrics).transpose()
            
            sample_metrics = pd.DataFrame(OrderedDict([
                                                       ('sample_id', [sample_id]), 
                                                       ('unpaired_mapped_reads', [df_metrics.ix[0, 'unpaired_reads_examined']]), 
                                                       ('paired_mapped_reads', [df_metrics.ix[0, 'read_pairs_examined']]), 
                                                       ('unpaired_duplicate_reads', [df_metrics.ix[0, 'unpaired_read_duplicates']]), 
                                                       ('paired_duplicate_reads', [df_metrics.ix[0, 'read_pair_duplicates']]), 
                                                       ('unmapped_reads', [df_metrics.ix[0, 'unmapped_reads']]), 
                                                       ('percent_duplicate_reads', [df_metrics.ix[0, 'percent_duplication']]), 
                                                       ('estimated_library_size', [df_metrics.ix[0, 'estimated_library_size']]), 
                                                       ]))
            
            metrics = metrics.append(sample_metrics)
    
    metrics = metrics.reset_index(drop=True)
    
    return metrics

def extract_alignment_metrics(dir):
    ''' Extract total reads with and without duplicates '''
    
    metrics_with_duplicates = pd.DataFrame()
    metrics_without_duplicates = pd.DataFrame()
    
    for file in os.listdir(dir):
        if file.endswith('.markdups.txt'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep='\t', header=6, skipfooter=2, engine='python')
            df.columns = [x.lower() for x in df.columns]
            
            total_reads = df.ix[2, 'total_reads']
            
            sample_metrics = pd.DataFrame({'sample_id': [sample_id], 'total_reads_with_duplicates': [total_reads]})
            
            metrics_with_duplicates = metrics_with_duplicates.append(sample_metrics)
            
        elif file.endswith('.rmdups.txt'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep='\t', header=6, skipfooter=2, engine='python')
            df.columns = [x.lower() for x in df.columns]
            
            total_reads = df.ix[2, 'total_reads']
            
            sample_metrics = pd.DataFrame({'sample_id': [sample_id], 'total_reads_without_duplicates': [total_reads]})
            
            metrics_without_duplicates = metrics_without_duplicates.append(sample_metrics)
    
    metrics = pd.merge(metrics_with_duplicates, metrics_without_duplicates, on='sample_id')
        
    return metrics

def extract_flagstat_metrics(dir):
    ''' Extract basic flagstat metrics '''
    
    metrics = pd.DataFrame()
    
    for file in os.listdir(dir):
        if file.endswith('.markdups.flagstat'):
            sample_id = file.split(".")[0]
            
            df = pd.read_table(os.path.join(dir, file), sep=r'\s\+\s0\s', header=None, names=['value', 'type'], engine='python')
            
            sample_metrics = pd.DataFrame(OrderedDict([
                                                       ('sample_id', [sample_id]), 
                                                       ('total_reads', [df.ix[0, 'value']]), 
                                                       ('total_mapped_reads', [df.ix[2, 'value']]), 
                                                       ('total_duplicate_reads', [df.ix[1, 'value']]), 
                                                       ('total_properly_paired', [df.ix[6, 'value']])
                                                       ]))
            
            metrics = metrics.append(sample_metrics)
    
    metrics = metrics.reset_index(drop=True)
    
    return metrics

#=======================================================================================================================
# Run script
#=======================================================================================================================

'''
NOTES: 
- All duplicate reads are mapped, since the duplicate flagging is done based on mapping coordinates!
- Alignment metrics not useful, as key information is contained in flagstat metrics

args.library_id = 'PX0218'
args.metrics_dir = '/share/lustre/asteif/projects/single_cell_indexing/alignment/PX0218/9_cycles/metrics'
args.out_file = '/share/lustre/asteif/projects/single_cell_indexing/test/PX0218.metrics_table.csv'

args.metrics_dir = '/share/scratch/asteif_temp/single_cell_indexing/merge/test_gatk_realign/metrics'
args.out_file = '/share/scratch/asteif_temp/single_cell_indexing/merge/test_gatk_realign/metrics/summary/test_gatk_realign.metrics_table.csv'

'''

def main():
    # alignment_metrics_dir = os.path.join(args.metrics_dir, 'alignment_metrics')
    duplication_metrics_dir = os.path.join(args.metrics_dir, 'duplication_metrics')
    flagstat_metrics_dir = os.path.join(args.metrics_dir, 'flagstat_metrics')
    insert_metrics_dir = os.path.join(args.metrics_dir, 'insert_metrics')
    wgs_metrics_dir = os.path.join(args.metrics_dir, 'wgs_metrics')
    
    # alignment_metrics = extract_alignment_metrics(alignment_metrics_dir)    
    duplication_metrics = extract_duplication_metrics(duplication_metrics_dir)
    flagstat_metrics = extract_flagstat_metrics(flagstat_metrics_dir)
    wgs_metrics = extract_wgs_metrics(wgs_metrics_dir)
    
    metrics_table = flagstat_metrics.merge(
                    duplication_metrics, on='sample_id').merge(
                    wgs_metrics, on='sample_id')
    
    if len(os.listdir(insert_metrics_dir)) > 0:
        insert_metrics = extract_insert_metrics(insert_metrics_dir)
        
        metrics_table = metrics_table.merge(insert_metrics, on='sample_id')
    
    if args.library_id:
        metrics_table.insert(0,'library_id', args.library_id)
    
    metrics_table.to_csv(args.out_file, index=False)

if __name__ == '__main__':
    main()
