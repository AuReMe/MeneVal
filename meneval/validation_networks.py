# -*- coding: utf-8 -*-
import os.path
import logging

from typing import List
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

def write_res(rxn, bact_pres, res_file: str):
    args_list = [rxn, bact_pres[0][0], round(bact_pres[0][1] * 100, 2), ';'.join(bact_pres[1])]
    args_list = [str(x) for x in args_list]
    with open(res_file, 'a') as f:
        f.write('\n' + '\t'.join(args_list))


# MAIN FUNCTION ========================================================================================================

def validation_networks(name_species: str, output: str, rxn_list: List[str], reactions_file: str,
                        group_file: str = None, group: str = None):
    init_logger(output)
    logging.info(f'Start searching for presence in {name_species} of reactions from Meneco output\n'
                 '======================================================================================\n')
    res_file = init_res_file(name_species, output)
    rxn_instance = create_rxn_instance(reactions_file, group_file, group)

    rxn_pres_dict = rxn_instance.get_rxn_presence(rxn_list)
    kept_rxn_set = set()
    for rxn, rxn_pres in rxn_pres_dict.items():
        if rxn_pres[0][0] != 0:
            write_res(rxn, rxn_pres, res_file)
            kept_rxn_set.add(rxn)

    logging.info(f'{len(kept_rxn_set)} reactions in {name_species} : {", ".join(kept_rxn_set)}\n'
                 f'Details in {res_file} file.')
    logging.info(f'\n================\n'
                 f'End of selection')
    return kept_rxn_set
