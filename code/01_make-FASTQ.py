import os
import logging
from datetime import (
  date,
  datetime
  )
import pandas as pd
from math import ceil
from shared import merge_tables
import subprocess

current_time = datetime.now()
date_and_time = current_time.strftime('%Y-%m-%d_%H-%M-%S')
date_and_time_pretty = current_time.strftime('%Y-%m-%d %H:%M:%S')

bcl_convert_version = get_software_version('bcl-convert')

def main():
    
    scop_id = snakemake.config['scop_id']
    com_id = snakemake.config['com_id']

    # path
    fastq_path = os.path.join(snakemake.config['project_path'],
                              scop_id,
                              snakemake.config['fastq_path'])
    
    bcl_data = merge_tables('reaction_sheet', 'reaction2library', 
                            'library_sheet',  'library2sequencing', 
                            'sequencing_sheet', snake_obj = snakemake.config)
    

    for bcl_folder, sub_sheet in bcl_data.groupby('bcl_folder'):  
        q_fastq_path = os.path.join(
          fastq_path,
          bcl_folder,
          bcl_convert_version)
        
        log_file = os.path.join(
          q_fastq_path,
          f'{date_and_time}.log')
        
        sheet_file = os.path.join(
          q_fastq_path,
          f'sample-sheet_{bcl_folder}.csv')
          
        os.makedirs(q_fastq_path, exist_ok=True)

        complete_file = os.path.join(
            q_fastq_path,
            'Logs',
            'FastqComplete.txt'
            )
            
        if os.path.isfile(complete_file):
          message = '################################################################################\n'
          message += f'FASTQ files already exists - skipping blc_folder: {bcl_folder}\n'
          message += '################################################################################\n'

          with open(os.path.realpath(os.path.join(q_fastq_path,f'_skipped_pipestance_{date_and_time}.log')), 'a') as log:
            log.writelines(message)
            print(message)
            continue

        logging.basicConfig(filename=log_file, level=logging.INFO,format='%(message)s')
        
        override_cycles = sub_sheet['override_cycles'].iloc[0]
        sequencing_id = sub_sheet['sequencing_id'].iloc[0]
                
        samplesheet = make_samplesheet(
            com_id = snakemake.config['com_id'], 
            indexes = sub_sheet['index'], 
            lanes = sub_sheet['lane'],
            bases_mask = override_cycles)
        
        
        with open(os.path.realpath(sheet_file), 'w') as sheet:
            sheet.writelines(samplesheet)
        
        bclconvert(
          scop_id = scop_id,
          com_id = com_id,
          sequencing_id = sequencing_id,
          bcl_folder = bcl_folder,
          override_cycles = override_cycles,
          fastq_path = fastq_path,
          bcl_data_path = snakemake.config['bcl_path'],
          threads = snakemake.threads,
          log_file = log_file)
          
    os.system(f'mkdir -p $(dirname {snakemake.output[0]}) && touch {snakemake.output[0]}')

def bclconvert(scop_id, com_id, sequencing_id, bcl_folder, override_cycles, fastq_path, bcl_data_path, threads, log_file):
    bcl_path = os.path.join(bcl_data_path, bcl_folder)
    fastq_folder = os.path.join(fastq_path, bcl_folder, bcl_convert_version)
    sample_sheet = os.path.join(fastq_path, bcl_folder, bcl_convert_version, f'sample-sheet_{bcl_folder}.csv')
    
    logging.info('#' * 80)
    logging.info('#####                                                                      #####')
    logging.info('#####                      Demultiplexing FASTQ files                      #####')
    logging.info('#####                                                                      #####')
    logging.info('#####             1/6                                                      #####')
    logging.info('#####           ▒▒▒▒▒▒▒▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░            #####')
    logging.info('#####                                                                      #####')
    logging.info('#' * 80)
    logging.info('#####')
    logging.info('##')
    logging.info(f'##\tSCOP ID:      \t\t{scop_id}')
    logging.info(f'##\tCOMUNEQAID ID:\t\t{com_id}')
    logging.info('##')
    logging.info(f'##\tsequencing_id:\t\t- {sequencing_id}')
    logging.info(f'##\tbcl_folder:   \t\t- {bcl_folder}')
    logging.info(f'##\toverride_cycles:\t- {override_cycles}')
    logging.info('##')
    logging.info('##')
    logging.info(f'##\t{date_and_time_pretty}')
    logging.info('##')
    logging.info('##')
    logging.info('###')
    logging.info('#' * 80)
    
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
