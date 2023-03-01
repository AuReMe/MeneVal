import os.path

import logging
import argparse
from typing import List, Set
from aucomana.utils.reactions import Reactions
from aucomana.utils.utils import get_grp_set


# CONSTANTS ============================================================================================================

# Command line arguments
def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, required=True, metavar='meneco output file', help='Meneco output file')
    parser.add_argument('-o', type=str, required=True, metavar='directory', help='The output directory')
    parser.add_argument('-rxn', type=str, required=True, metavar='reactions.tsv file',
                        help='reactions.tsv file from compare-padmets output')
    parser.add_argument('-n', type=str, required=True, metavar='source name', help='Name of source')
    parser.add_argument('-gf', type=str, required=False, metavar='groups template file',
                        help='The groups template tsv file from AuCoMe runs')
    parser.add_argument('-g', type=str, required=False, metavar='name of group', help='The name of the group to select')
    args = parser.parse_args()
    if args.g is None:
        args.g = ''
    return args.m, args.o, args.rxn, args.n, args.gf, args.g


MENECO_TSV, OUTPUT, RXN, NAME, GRP_FILE, GRP = get_command_line_args()


# CREATE REACTION INSTANCE

if GRP_FILE is not None:
    sp_l = list(get_grp_set(group_file=GRP_FILE, group=GRP))
    RXN = Reactions(file_reactions_tsv=RXN, species_list=sp_l)
else:
    RXN = Reactions(file_reactions_tsv=RXN)


# RES
RES_FILE = os.path.join(OUTPUT, 'res_val2.tsv')
RES_HEADER = ['RXN', f'Nb {NAME} presence', f'{NAME} presence %', f'{NAME} list']
with open(RES_FILE, 'w') as F:
    F.write('\t'.join(RES_HEADER))

# Init logger
LOG_FILE = os.path.join(OUTPUT, 'meneco_validation.log')
logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(message)s')


# FUNCTIONS ============================================================================================================

def extract_rxn_from_meneco() -> List[str]:
    """Extract the list of reactions corresponding to the union of reactions from solutions found by Meneco.
    The extraction is done from a meneco_output.tsv file created with padmet enhance_meneco_output.

    Returns
    -------
    List[str]
        List of reactions
    """
    logging.info('Extracting reactions from Meneco output')
    rxn_list = list()
    with open(MENECO_TSV, 'r') as meneco:
        for line in meneco:
            rxn = line.split('\t')[0]
            rxn_list.append(rxn)
    del rxn_list[0]
    logging.info(f'Total of {len(rxn_list)} reactions\n')
    return rxn_list


def write_res(rxn, bact_pres):
    args_list = [rxn, bact_pres[0][0], round(bact_pres[0][1] * 100, 2), ';'.join(bact_pres[1])]
    args_list = [str(x) for x in args_list]
    with open(RES_FILE, 'a') as f:
        f.write('\n' + '\t'.join(args_list))


def create_new_meneco_tsv(kept_rxn: Set[str]):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    new_meneco = os.path.join(OUTPUT, 'meneco_output_filtered.tsv')
    with open(MENECO_TSV, 'r') as in_meneco, open(new_meneco, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                line[-2] = f'Potential {NAME} source'
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {new_meneco}')


# MAIN FUNCTION ========================================================================================================

def meneco_validation_2():
    logging.info(f'Start searching for presence in {NAME} of reactions from Meneco output\n'
                 '======================================================================================\n')

    rxn_l = extract_rxn_from_meneco()
    rxn_pres_dict = RXN.get_rxn_presence(rxn_l)
    kept_rxn_set = set()
    for rxn, rxn_pres in rxn_pres_dict.items():
        if rxn_pres[0][0] != 0:
            write_res(rxn, rxn_pres)
            kept_rxn_set.add(rxn)

    logging.info(f'{len(kept_rxn_set)} reactions in {NAME} : {", ".join(kept_rxn_set)}\n'
                 f'Details in {RES_FILE} file.')
    create_new_meneco_tsv(kept_rxn_set)
    logging.info(f'\n================\n'
                 f'End of selection')


meneco_validation_2()
