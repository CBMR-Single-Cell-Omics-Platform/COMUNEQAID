import os
import logging
from datetime import date
import pandas as pd
from math import ceil
from shared import merge_tables

def main():
    # Initialize logging
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0], 
                        format='%(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

    print('#' * 80)
    print('#####                                                                      #####')
    print('#####                      Demultiplexing FASTQ files                      #####')
    print('#####                                                       ** / *****     #####')
    print('#' * 80)
    print('#####\n##')
    print('##\tDemultiplexing libraries:')
    
    fastq_path = os.path.join(snakemake.config['project_path'],
                              snakemake.config['scop_id'],
                              snakemake.config['fastq_path'])
    
    os.makedirs(fastq_path, exist_ok=True)
    
    
    bcl_data = merge_tables('reaction_sheet', 'reaction2library', 
                            'library_sheet',  'library2sequencing', 
                            'sequencing_sheet', snake_obj = snakemake.config)

    for bcl_folder, sub_sheet in bcl_data.groupby('bcl_folder'):  
        sheet_file = os.path.join(fastq_path,
                                  f'{bcl_folder}_SampleSheet.csv')
        
        log_path = os.path.join(snakemake.config['project_path'],
                                snakemake.config['scop_id'],
                                snakemake.config['log_path'])
        os.makedirs(log_path, exist_ok=True)
        
        log_file = os.path.join(log_path, \
                                f'{bcl_folder}_bcl-convert.log')
        
        complete_file = os.path.join(
            fastq_path,
            bcl_folder,
            'Logs',
            'FastqComplete.txt'
            )
        if os.path.isfile(complete_file):
            print('##\t-\tpipestance exists - skipping')
            continue
        
        override_cycles = sub_sheet['override_cycles'].iloc[0]
                
        samplesheet = make_samplesheet(
            com_id = snakemake.config['com_id'], 
            indexes = sub_sheet['index'], 
            lanes = sub_sheet['lane'],
            bases_mask = override_cycles)
        
        
        with open(os.path.realpath(sheet_file), 'w') as sheet:
            sheet.writelines(samplesheet)
        
        
        bclconvert(bcl_folder = bcl_folder,
                   fastq_path = fastq_path,
                   bcl_data_path = snakemake.config['bcl_path'],
                   threads = snakemake.threads,
                   log_file = log_file)
    os.system(f'mkdir -p $(dirname {snakemake.output[0]}) && touch {snakemake.output[0]}')

def bclconvert(bcl_folder, fastq_path, bcl_data_path, threads, log_file):
    bcl_path = os.path.join(bcl_data_path, bcl_folder)
    fastq_folder = os.path.join(fastq_path, bcl_folder)
    sample_sheet = os.path.join(fastq_path, f'{bcl_folder}_SampleSheet.csv')
    
    os.system(
        f'bcl-convert '
        f'--output-directory {fastq_folder} '
        f'--sample-sheet {sample_sheet} '
        f'--bcl-num-conversion-threads {threads} '
        f'--bcl-num-compression-threads {threads} '
        f'--bcl-num-decompression-threads {ceil(threads * 0.5)} '
        f'--bcl-input-directory {bcl_path} '
        f'--force '
        f'2>&1 | tee {log_file} -a'
    )


def make_samplesheet(com_id, indexes, lanes, bases_mask = None):
    override = '' if (bases_mask == '') else f'OverrideCycles,{bases_mask}\n'
    mergelanes = 'TRUE' if all(lanes == '*') else 'FALSE'
    samplesheet = (
        f'[Header]\n'
        f'Date,{date.today()}\n'
        f'FileFormatVersion,2\n'
        f'Project Name,{com_id}\n'
        f'\n'
        f'[BCLConvert_Settings]\n'
        f'NoLaneSplitting,{mergelanes}\n'
        f'{override}'
        f'\n'
        f'[BCLConvert_Data]\n'
        )
    if (all(lanes == '*')):
        samplesheet += 'Sample_ID,index,index2\n'
        for idx in indexes:
            samplesheet += (f'{idx},{all_indexes[idx][0]},'
                            f'{all_indexes[idx][1]}\n')
    else:
        samplesheet += 'Lane,Sample_ID,index,index2\n'
        for idx, lane in zip(indexes, lanes):
            samplesheet += (f'{lane},{idx},{all_indexes[idx][0]},'
                            f'{all_indexes[idx][1]}\n')
    return samplesheet

def index_reader(file, idx2 = 2, skip = None):
    idx_csv = pd.read_csv(file, skiprows = skip, usecols=[0,1,idx2])
    idx_csv.columns = ['idx', 'seq1', 'seq2']
    indexes = {row.idx : [row.seq1, row.seq2] for _, row in idx_csv.iterrows()}
    return indexes

all_indexes = index_reader('data/Dual_Index_Kit_TT_Set_A.csv', 3, 3) | \
    index_reader('data/TruSeq_I7_indexes.csv')


main()
