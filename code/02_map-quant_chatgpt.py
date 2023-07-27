import os
import glob
import logging
from collections import defaultdict
import shared


def main():
    # Initialize logging
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0], format='%(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('#' * 80)
    logging.info('#####                                                                      #####')
    logging.info('#####                      Mapping and quantification                      #####')
    logging.info('#####                                                      *** / *****     #####')
    logging.info('#' * 80)
    logging.info('#####\n##')
    logging.info('##\tMapping and quantifying 10x libraries:')

    indexes_10x = index_list(snakemake.config['pin_10x'])
    for index, seq_names in indexes_10x.items():
        process_index(index, seq_names, '10x')

    if snakemake.config['workflow'] == '10x + HTO':
        indexes_hto = index_list(snakemake.config['pin_hto'])
        for index, seq_names in indexes_hto.items():
            process_index(index, seq_names, 'hto')

    os.system(f'touch {snakemake.output}')


def index_list(pin_dict):
    indexes = defaultdict(list)
    for seq_name in pin_dict:
        for index in pin_dict[seq_name]:
            indexes[index].append(seq_name)
    return indexes


def prepare_hto_index(index, path_quant_out):
    features_tsv = os.path.join(path_quant_out, 'features.tsv')
    t2g = os.path.join(path_quant_out, 't2g_hto.tsv')

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


def salmon_alevin(files_read_1, files_read_2, index_path, salmon_alevin_param, path_quant_out):
    read_1_string = ' '.join(files_read_1)
    read_2_string = ' '.join(files_read_2)
    logging.info('##\tmapping with alevin')
    os.system(
        'salmon alevin '
        '-l ISR '
        f'-i {index_path} '
        f'-1 {read_1_string} '
        f'-2 {read_2_string} '
        f'-p {snakemake.threads} '
        f'-o {path_quant_out}/map '
        '--sketch '
        f'{salmon_alevin_param} '
        f'2>&1 | tee {path_quant_out}/alevin-fry.log -a'
    )


def alevin_fry(path_quant_out, direction, tgmap):
    logging.info('##\tgenerating permit list')
    whitelist = snakemake.config['whitelist']
    os.system(
        'alevin-fry generate-permit-list '
        f'-d {direction} '
        f'-i {path_quant_out}/map '
        f'-o {path_quant_out}/quant '
        f'-u {whitelist} '
        f'2>&1 | tee {path_quant_out}/alevin-fry.log -a'
    )

    logging.info('##\tcollating')
    os.system(
        'alevin-fry collate '
        f'-r {path_quant_out}/map '
        f'-i {path_quant_out}/quant '
        '-t 1 '
        f'2>&1 | tee {path_quant_out}/alevin-fry.log -a'
    )
    logging.info('##\tquantifying')
    os.system(
        'alevin-fry quant '
        f'--threads {snakemake.threads} '
        f'--input-dir {path_quant_out}/quant '
        f'--output-dir {path_quant_out}/res '
        '--resolution cr-like '
        '--use-mtx '
        f'--tg-map {tgmap} '
        f'2>&1 | tee {path_quant_out}/alevin-fry.log -a'
    )


def count_lines(file):
    with open(file) as f:
        return sum(1 for _ in f)


def process_index(index, seq_names, index_type):
    logging.info(f'##\t-\tprocessing index:\t{index}')
    path_quant_out = os.path.join(
        snakemake.config['out_path'],
        snakemake.config['com_id'],
        index_type,
        index
    )

    seq_names_join = ', '.join(seq_names)
    logging.info(f'##\t\t(from {seq_names_join})')

    run_complete = os.path.isdir(os.path.join(path_quant_out, "res"))

    if run_complete:
        logging.info('##\t-\tpipestance exists - skipping')
        return

    os.makedirs(path_quant_out, exist_ok=True)

    files_read_1 = []
    files_read_2 = []
    for seq_name in seq_names:
        flowID = seq_name.split('_')[3][1:]
        fastq_folder = os.path.join(snakemake.config['fastq_path'], seq_name, 'fastq')
        files_read_1 += glob.glob(f'{fastq_folder}/?{flowID}/*{index}*R1*.fastq.gz')
        files_read_2 += glob.glob(f'{fastq_folder}/?{flowID}/*{index}*R2*.fastq.gz')

    if index_type == 'hto':
        index_path, alevin_fry_tgmap = prepare_hto_index(index, path_quant_out)
        direction = 'fw'
        salmon_alevin_param = (
            '--read-geometry 2[1-15] '
            '--bc-geometry 1[1-16] '
            '--umi-geometry 1[17-26]'
        )
        usa_mode = 'false'
    else:
        index_path = snakemake.config['salmon_index']
        direction = 'both'
        alevin_fry_tgmap = snakemake.config['T3G']
        t2g = snakemake.config['T3G']
        salmon_alevin_param = f'--chromiumV3 --tgMap {t2g}'
        usa_mode = 'true'

    salmon_alevin(files_read_1, files_read_2, index_path, salmon_alevin_param, path_quant_out)
    alevin_fry(path_quant_out, direction, alevin_fry_tgmap)

    alevin_res = os.path.join(path_quant_out, 'res', 'alevin')
    count_genes = count_lines(f'{alevin_res}/quants_mat_cols.txt')
    count_cells = count_lines(f'{alevin_res}/quants_mat_rows.txt')

    meta_info = {
        'alt_resolved_cell_numbers': [],
        'dump_eq': False,
        'num_genes': count_genes,
        'num_quantified_cells': count_cells,
        'resolution_strategy': 'CellRangerLike',
        'usa_mode': usa_mode
    }

    with open(f'{path_quant_out}/res/meta_info.json', 'w') as f:
        json.dump(meta_info, f)

hto_sequences = pd.read_csv('data/totalseq-a-hashtags.csv',
                            delimiter='\t',
                            usecols=['Barcode', 'Sequence'])

main()


pd.merge(sample_data, hto_sequences, left_on='hto', right_on='Barcode')

sample_data['hto'].isin(hto_sequences['Barcode'])
