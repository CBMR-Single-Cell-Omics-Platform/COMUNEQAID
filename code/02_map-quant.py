# Imports
import os
import glob
import logging
import pandas as pd
from shared import merge_tables

def main():
    # Init
    logging.basicConfig(level=logging.INFO,
                    filename=snakemake.log[0],
                    format='%(message)s')
                    
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('#'*80)
    logging.info('#####                                                                      #####')
    logging.info('#####                      Mapping and quantification                      #####')
    logging.info('#####                                                      *** / *****     #####')
    logging.info('#'*80)
    logging.info('#####\n##')
    logging.info('##\tMapping and quantifying 10x libraries:')
    
    fastq_path = os.path.join(snakemake.config['project_path'],
                              snakemake.config['scop_id'],
                              snakemake.config['fastq_path'])
    
    quant_path = os.path.join(snakemake.config['project_path'],
                            snakemake.config['scop_id'],
                            snakemake.config['out_path'],
                            snakemake.config['com_id']
                            )
    
    sample_data = merge_tables(
        'sample_sheet', 'sample2reaction', 
        'reaction_sheet', snake_obj = snakemake.config
        )
    
    reaction_data = merge_tables(
        'reaction_sheet', 'reaction2library', 
        'library_sheet',  'library2sequencing', 
        'sequencing_sheet', snake_obj = snakemake.config
        )
        
    # Prepare output directories
    for reaction in set(reaction_data['reaction_id']):
        reaction_path = os.path.join(quant_path, reaction)
        os.makedirs(reaction_path, exist_ok=True)
    
        
    # Align
    for idx, reaction in reaction_data.groupby(['reaction_id', 'library_type']):
        logging.info(f'##\t-\tproccessing {idx[1]} for reaction:\t{idx[0]}')
        
        out_folder = os.path.join(
            quant_path, idx[0], idx[1]
            )
        os.makedirs(out_folder, exist_ok=True)
        
        meta_info_file = os.path.join(
            out_folder, 'res', 'meta_info.json'
            )
        
        if os.path.isfile(meta_info_file):
            print('##\t-\tpipestance exists - skipping')
            continue
        
        files_read_1 = list()
        files_read_2 = list()
        
        def helper(x, row):
            lane_string = '' if row['lane'] == '*' else f'_L00{row["lane"]}'
            file = row['index'] + '_S*' + lane_string + x
            return os.path.join(fastq_path, row['bcl_folder'], file)
        
        for _, row in reaction.iterrows():
            files_read_1 += glob.glob(helper('_R1_*.fastq.gz', row))
            files_read_2 += glob.glob(helper('_R2_*.fastq.gz', row))
        
        library_type = reaction['library_type'].iloc[0]
        
        if library_type == '10x':
            align = reaction['align'].iloc[0]
            index_path = snakemake.config['salmon_index'][align]
            t3g = snakemake.config['T3G'][align]
            t2g = snakemake.config['T2G'][align]
        elif library_type == 'hto':
            samples = sample_data.loc[sample_data['reaction_id'] == idx[0]]
            hash_codes = pd.merge(samples, hto_sequences, 
                                 left_on='hto', right_on='Barcode')
            hash_codes = hash_codes[['sample_id', 'Sequence']]
            index_path, t3g = prepare_hto_index(hash_codes, out_folder)
            t2g = ''
        else:
            raise Exception(f"Library type {library_type} not supported")
        
        process_index(  
            files_read_1 = files_read_1, 
            files_read_2 = files_read_2, 
            library_type = library_type,
            index_path = index_path,
            t2g = t2g,
            alevin_fry_tgmap = t3g,
            path_quant_out = out_folder,
            meta_info_file = meta_info_file)
        
    os.system(f'mkdir -p $(dirname {snakemake.output}) && touch {snakemake.output}')
    

def process_index(files_read_1, files_read_2, library_type, index_path, 
                  t2g, alevin_fry_tgmap, path_quant_out, meta_info_file):
    
    if library_type == "hto":
        direction = "fw"
        salmon_alevin_param = f'\
            --read-geometry 2[1-15] \
            --bc-geometry 1[1-16] \
            --umi-geometry 1[17-26]'
        usa_mode = 'false'
    else:
        direction = "both"
        salmon_alevin_param = f'\
            --chromiumV3 \
            --tgMap {t2g}'
        usa_mode = 'true'

    salmon_alevin(files_read_1, files_read_2, index_path, 
                  salmon_alevin_param, path_quant_out)
    alevin_fry(path_quant_out, direction, alevin_fry_tgmap)
    
    alevin_res = os.path.join(path_quant_out, 'res', 'alevin')
    count_genes = count_lines(os.path.join(alevin_res, 'quants_mat_cols.txt'))
    count_cells = count_lines(os.path.join(alevin_res, 'quants_mat_rows.txt'))

    with open(meta_info_file, 'w') as f:
            f.write('{')
            f.write('  "alt_resolved_cell_numbers": [],')
            f.write('  "dump_eq": false,')
            f.write(f'  "num_genes": {count_genes},')
            f.write(f'  "num_quantified_cells": {count_cells},')
            f.write('  "resolution_strategy": "CellRangerLike",')
            f.write(f'  "usa_mode": {usa_mode}')
            f.write('}')


def prepare_hto_index(hash_codes, path_quant_out):
    features_tsv = os.path.join(path_quant_out, 'features.tsv')
    t2g = os.path.join(path_quant_out, 't2g_hto.tsv')
    
    with open(features_tsv, 'w') as features, open(t2g, 'w') as t2g_file:
        for _, row in hash_codes.iterrows():
            print(row['sample_id'] + '\t' + row['Sequence'], file=features)
            print(row['sample_id'] + '\t' + row['sample_id'], file=t2g_file)
    

    # The transcript2gene file for the HTO is a 1:1 mapping between the feature names
    with open(features_tsv, 'r') as in_file, open(t2g, 'w') as out_file:
        for line in in_file:
            feature_id = line.split('\t')[0]
            print(feature_id + '\t' + feature_id, file=out_file)

    logging.info('##\tbuilding HTO index')
    index_path = os.path.join(path_quant_out, 'hto_index')
    os.system(
        'salmon index '
        f'-t {features_tsv} '
        f'-i {index_path} '
        '--features '
        '-k 7 '
        f'2>&1 | tee {path_quant_out}/alevin-fry.log -a'
    )
    return index_path, t2g

def salmon_alevin(files_read_1, files_read_2, index_path, 
                  salmon_alevin_param, path_quant_out):
    read_1_string = ' '.join(files_read_1)
    read_2_string = ' '.join(files_read_2)
    logging.info(f'##\tmapping with alevin')
    os.system(f'salmon alevin \
        -l ISR \
        -i {index_path} \
        -1 {read_1_string} \
        -2 {read_2_string} \
        -p {snakemake.threads} \
        -o {path_quant_out}/map \
        --sketch \
        {salmon_alevin_param} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')

def alevin_fry(path_quant_out, direction, tgmap):
    logging.info(f'##\tgenerating permit list')
    whitelist = snakemake.config['whitelist']
    os.system(f'alevin-fry generate-permit-list \
        -d {direction} \
        -i {path_quant_out}/map \
        -o {path_quant_out}/quant \
        -u {whitelist} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    
    logging.info(f'##\tcollating')
    os.system(f'alevin-fry collate \
        -r {path_quant_out}/map \
        -i {path_quant_out}/quant \
        -t 1 \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')
    logging.info(f'##\tquantifying')
    os.system(f'alevin-fry quant \
        --threads {snakemake.threads} \
        --input-dir {path_quant_out}/quant \
        --output-dir {path_quant_out}/res \
        --resolution cr-like \
        --use-mtx \
        --tg-map {tgmap} \
        2>&1 | tee {path_quant_out}/alevin-fry.log -a')

def count_lines(file):
    i = -1
    with open(file) as f:
            for i, _ in enumerate(f):
                pass
    return i + 1

hto_sequences = pd.read_csv('data/totalseq-a-hashtags.csv',
                            delimiter='\t',
                            usecols=['Barcode', 'Sequence'])



main()
