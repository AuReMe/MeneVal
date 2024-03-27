import os.path
import logging

from typing import List, Tuple, Set

import aucomana.utils.reactions
from aucomana.utils.reactions import Reactions
from aucomana.utils.utils import get_grp_set


# ENVIRONMENT ==========================================================================================================

# RES
def init_res_file(name_species: str, output: str) -> str:
    """ Initialize the results file and return its path

    Parameters
    ----------
    name_species: str
        Name of the species studied (class of species per example)
    output: str
        Output directory to store the results file

    Returns
    -------
    str
        Path to the results file
    """
    res_file = os.path.join(output, f'{name_species}_res_validation_networks.tsv')
    res_header = ['RXN', f'Nb {name_species} presence', f'{name_species} presence %', f'{name_species} list']
    with open(res_file, 'w') as F:
        F.write('\t'.join(res_header))
    return res_file


# Init logger
def init_logger(output: str):
    """ Initialize logger file

    Parameters
    ----------
    output: str
        Path to the output directory to store the log file
    """
    log_file = os.path.join(output, 'networks_validation.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s', force=True)


# FUNCTIONS ============================================================================================================

def create_rxn_instance(reactions_file: str, group_file: str or None, group: str or None) -> \
        aucomana.utils.reactions.Reactions:
    """ Create an aucomana Reactions instance from reactions.tsv file created from comparison of padmet networks.
    Selects a specified group of species if group_file and group parameters filled.

    Parameters
    ----------
    reactions_file: str
        Path to reactions.tsv file created from comparison of padmet networks
    group_file: str or None
        Path to group_template.tsv file from aucome analysis step, if None all species will be selected
    group: str or None
        Name of the group to select, if None all species will be selected

    Returns
    -------
    aucomana.utils.reactions.Reactions
        Reactions instance from aucomana package
    """
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


def write_res(rxn: str, rxn_pres: Tuple[Tuple[int, float], Set[str]], res_file: str):
    """ Writes the line of the result for a given reaction in the result file.

    Parameters
    ----------
    rxn: str
        Reaction to consider
    rxn_pres: tuple[tuple[int, float], set[str]]
        Indication of the presence of the reaction among the species :
        tuple[tuple[number_of_species_having_the_reaction, percentage_of_species_having_the_reaction],
              list[species_having_the_reaction]]
    res_file: str
        Path to the result file to write the result
    """
    args_list = [rxn, rxn_pres[0][0], round(rxn_pres[0][1] * 100, 2), ';'.join(rxn_pres[1])]
    args_list = [str(x) for x in args_list]
    with open(res_file, 'a') as f:
        f.write('\n' + '\t'.join(args_list))


# MAIN FUNCTION ========================================================================================================

def validation_networks(name_species: str, output: str, rxn_list: List[str], reactions_file: str,
                        group_file: str = None, group: str = None) -> Set[str]:
    """ Checks from a list of reactions if some are present in networks of a group of species. The check is done from
    a comparison file of padmet networks.

    Parameters
    ----------
    name_species: str
        Name of the species studied (class of species per example)
    output: str
        Output directory to store the results files
    rxn_list: list[str]
        List of reactions to check
    reactions_file: str
        Path to reactions.tsv file created from comparison of padmet networks
    group_file: str or None
        Path to group_template.tsv file from aucome analysis step, if None all species will be selected
    group: str or None
        Name of the group to select, if None all species will be selected

    Returns
    -------
    set[str]
        Set of reactions that were found in other group species networks
    """
    init_logger(output)
    logging.info(f'Start searching for presence in {name_species} of reactions'
                 '==================================================================\n')
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
