# -*- coding: utf-8 -*-
import os.path
import logging

from typing import List, Set
from aucomana.utils.reactions import Reactions
from aucomana.utils.utils import get_grp_set


# ENVIRONMENT ==========================================================================================================

def create_rxn_instance(reactions_file: str, group_file: str, group: str):
    if group_file is not None:
        if group is None:
            logging.warning('Group file specified but no group name was given, will consider all species.')
            rxn = Reactions(file_reactions_tsv=reactions_file)
            return rxn
        else:
            sp_l = list(get_grp_set(group_file=group_file, group=group))
            rxn = Reactions(file_reactions_tsv=reactions_file, species_list=sp_l)
            return rxn
    rxn = Reactions(file_reactions_tsv=reactions_file)
    return rxn


# RES
def init_res_file(name_species: str, output: str):
    res_file = os.path.join(output, 'res_val2.tsv')
    res_header = ['RXN', f'Nb {name_species} presence', f'{name_species} presence %', f'{name_species} list']
    with open(res_file, 'w') as F:
        F.write('\t'.join(res_header))
    return res_file


# Init logger
def init_logger(output: str):
    log_file = os.path.join(output, 'networks_validation.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')


# FUNCTIONS ============================================================================================================

def extract_rxn_from_meneco(meneco_tsv: str) -> List[str]:
    """Extract the list of reactions corresponding to the union of reactions from solutions found by Meneco.
    The extraction is done from a meneco_output.tsv file created with padmet enhance_meneco_output.

    Returns
    -------
    List[str]
        List of reactions
    """
    logging.info('Extracting reactions from Meneco output')
    rxn_list = list()
    with open(meneco_tsv, 'r') as meneco:
        for line in meneco:
            rxn = line.split('\t')[0]
            rxn_list.append(rxn)
    del rxn_list[0]
    logging.info(f'Total of {len(rxn_list)} reactions\n')
    return rxn_list


def write_res(rxn, bact_pres, res_file: str):
    args_list = [rxn, bact_pres[0][0], round(bact_pres[0][1] * 100, 2), ';'.join(bact_pres[1])]
    args_list = [str(x) for x in args_list]
    with open(res_file, 'a') as f:
        f.write('\n' + '\t'.join(args_list))


def create_new_meneco_tsv(kept_rxn: Set[str], name_species: str, output: str, meneco_tsv: str):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    name_species : str
        type of species considered
    output : str
        path to the output folder
    meneco_tsv : str
        Meneco results in TSV format (output from enhanced_meneco_output in padmet)
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    new_meneco = os.path.join(output, 'meneco_output_filtered.tsv')
    with open(meneco_tsv, 'r') as in_meneco, open(new_meneco, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                line[-2] = f'Potential {name_species} source'
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {new_meneco}')


# MAIN FUNCTION ========================================================================================================

def validation_networks(name_species: str, output: str, meneco_tsv: str, reactions_file: str,
                               group_file: str = None, group: str = None):
    init_logger(output)
    logging.info(f'Start searching for presence in {name_species} of reactions from Meneco output\n'
                 '======================================================================================\n')
    res_file = init_res_file(name_species, output)
    rxn_instance = create_rxn_instance(reactions_file, group_file, group)

    rxn_l = extract_rxn_from_meneco(meneco_tsv)
    rxn_pres_dict = rxn_instance.get_rxn_presence(rxn_l)
    kept_rxn_set = set()
    for rxn, rxn_pres in rxn_pres_dict.items():
        if rxn_pres[0][0] != 0:
            write_res(rxn, rxn_pres, res_file)
            kept_rxn_set.add(rxn)

    logging.info(f'{len(kept_rxn_set)} reactions in {name_species} : {", ".join(kept_rxn_set)}\n'
                 f'Details in {res_file} file.')
    create_new_meneco_tsv(kept_rxn_set, name_species, output, meneco_tsv)
    logging.info(f'\n================\n'
                 f'End of selection')
