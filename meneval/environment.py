"""
Functions to initialize and check the input files and folder required for meneco_validation workflow.
The initialisation will create all required folders.
The check will verify if all the required files are placed in the corrects folders.
"""
import os
import logging
import shutil
from typing import Dict, List, Set

# BASE ENVIRONMENT =====================================================================================================

# MAIN DIRECTORIES
INPUT = os.path.join('Input')
OUTPUT = os.path.join('Output')

ENRICH_D = 'Enrichment'
DATABASE_D = 'DataBase'
NETWORK_D = 'Networks'
SEEDS_D = 'Seeds'
SPECIES_D = 'Species_seq'
TARGETS_D = 'Targets'
BLASTP_D = 'BlastP'
MENECO_D = 'Meneco'

# TYPE DIRECTORIES
TSV_D = 'TSV'
SBML_D = 'SBML'
PADMET_D = 'PADMET'

FILTERED_D = 'Filtered_TSV'
TOOL_OUTPUTS_D = 'Json_outputs'

# FILES
REACTIONS_TSV = 'reactions.tsv'
GROUPS_TSV = 'group_template.tsv'
IN_TARGETS = 'targets.tsv'
IN_SEEDS = 'seeds.tsv'
IN_ARTEFACTS = 'artefacts.tsv'

# EXTENSIONS
TSV_EXT = '.tsv'
PADMET_EXT = '.padmet'
FASTA_EXT = '.fasta'
FAA_EXT = '.faa'
SBML_EXT = '.sbml'

# STEPS
BLASTP = 'BLASTP'
ENRICH = 'ENRICHMENT'
FILL = 'FILL'
EXCLUDE_E = 'NO_ENRICH_RXN'

# STR
GROUP_ALL = 'ALL'

# UTILITY FUNCTIONS ====================================================================================================
def get_file_from_ext(path: str, ext: str):
    """ Returns the path of the unique file with the given extension in the given path.

    Parameters
    ----------
    path : str
        The path where to find the file with the extension
    ext : str
        The extension of the file to found in the path

    Returns
    -------
    str
        The path of the unique file with the given extension in the given path

    Raises
    ------
    OSError
        If 0 or +1 files found with the given extension in the given path
    """
    files = [x for x in os.listdir(path) if x.endswith(ext)]
    if len(files) == 0:
        logging.info(f'No file with the extension {ext} in path {path}')
    elif len(files) > 1:
        logging.info(f'More than 1 file with the extension {ext} in path {path}')
    else:
        return os.path.join(path, files[0])


# INPUT FILES TO CREATE ================================================================================================

# SEEDS - ARTEFACTS

SEEDS_TSV = os.path.join(INPUT, SEEDS_D, f'seeds_medium{TSV_EXT}')
SEEDS_ARTEFACTS = {TSV_D: os.path.join(INPUT, SEEDS_D, f'seeds_artefacts{TSV_EXT}'),
                   SBML_D: os.path.join(INPUT, SEEDS_D, f'seeds_artefacts{SBML_EXT}')}

# TARGETS - BIOMASS

BIOMASS_TSV = os.path.join(INPUT, TARGETS_D, f'biomass{TSV_EXT}')
TARGETS = {TSV_D: os.path.join(INPUT, TARGETS_D, f'temp_targets{TSV_EXT}'),
           SBML_D: os.path.join(INPUT, TARGETS_D, f'targets{SBML_EXT}')}

# DATABASE

DB_SBML = os.path.join(INPUT, DATABASE_D, f'database{SBML_EXT}')

# NETWORKS

# START NETWORKS

MEDIUM_NW = os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'0_medium{PADMET_EXT}')
BASE_NW = {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'0_base{PADMET_EXT}'),
           SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'0_base{SBML_EXT}')}


def get_file_comp(file: str) -> List:
    """ Returns the decomposition of elements of a file name to the list of elements
    Examples: 1_Group1_ENRICH.padmet -> [1, Group1, ENRICH, padmet]

    Parameters
    ----------
    file: str
        Name of the file

    Returns
    -------
    List
        List of file name elements
    """
    file = file.split('.')
    name = file[0].split('_')
    name[0] = int(name[0])
    return name + [file[1]]


def get_num(step: str, group: str = None) -> int:
    """ Return the number of the step basing on number of step already done. If the current step files are found,
    return the number given in the already done step.

    Parameters
    ----------
    step: str
        Name of the step to find number associated with
    group: str (optional)
        Group name for enrichment step

    Returns
    -------
    int
        Number of the step.
    """
    num = len(os.listdir(os.path.join(OUTPUT, NETWORK_D, PADMET_D)))
    for file in os.listdir(os.path.join(OUTPUT, NETWORK_D, PADMET_D)):
        file_comp = get_file_comp(file)
        if step == ENRICH and file_comp[1] == group:
            return file_comp[0]
        elif step != ENRICH and file_comp[-2] == step:
            return file_comp[0]
    return num


# GAPFILLING NETWORKS
def get_nw_path(step, group=None) -> Dict[str, str]:
    """ Returns Network files name for PADMET and SBML extensions in a dictionary.

    Returns
    -------
    Dict[str, str]
        Network files name for PADMET and SBML extensions.
    """
    num = get_num(step, group)
    if step == ENRICH:
        return {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'{num}_{group}_{ENRICH}{PADMET_EXT}'),
                SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'{num}_{group}_{ENRICH}{SBML_EXT}')}
    else:
        return {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'{num}_{step}{PADMET_EXT}'),
                SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'{num}_{step}{SBML_EXT}')}


# INITIALIZATION =======================================================================================================

def create_folders():
    """ Creates all the input and outputs folders required for the workflow
    """
    logging.info('Running init step : creating directories :\n'
                 '=========================================\n')

    dir_archi = {INPUT: [DATABASE_D,
                         ENRICH_D,
                         SPECIES_D,
                         NETWORK_D,
                         SEEDS_D,
                         TARGETS_D,
                         ],

                 OUTPUT: [BLASTP_D,
                          ENRICH_D,
                          {NETWORK_D: [PADMET_D,
                                       SBML_D]},
                          {MENECO_D: [FILTERED_D,
                                      TOOL_OUTPUTS_D,
                                      TSV_D]}
                          ]
                 }
    create_dir_rec(dir_archi)

    logging.info('\n--------------\nInit step done\n')


def create_dir_rec(dir_dict: Dict[str, List[str or Dict[...]]], path: str = ''):
    """ Create directories from dictionary in a path recursively

    Parameters
    ----------
    dir_dict : Dict
        The dictionary of folders architecture
    path : str optional, default=''
         The parent path where to create folders
    """
    for parent_rep, child_list in dir_dict.items():
        parent_path = os.path.join(path, parent_rep)
        if not os.path.exists(parent_path):
            os.mkdir(parent_path)
            logging.info(f'{parent_path} created')
        else:
            logging.info(f'{parent_path} already exists')
        for child in child_list:
            if type(child) == str:
                child_path = os.path.join(path, parent_rep, child)
                if not os.path.exists(child_path):
                    os.mkdir(child_path)
                    logging.info(f'{child_path} created')
                else:
                    logging.info(f'{child_path} already exists')
            else:
                create_dir_rec(child, os.path.join(path, parent_rep))


def get_enrich_groups() -> Set[str]:
    """ Returns the set of groups for enrichment step according to folders name in the input enrichment folder.

    Returns
    -------
    Set[str]
        Set of groups for enrichment step
    """
    enrich_input_dir = os.path.join(INPUT, ENRICH_D)
    groups = set()
    if os.listdir(enrich_input_dir) != list():
        for directory in os.listdir(enrich_input_dir):
            groups.add(directory)
    return groups


def get_enrich_reactions_files() -> Dict[str, str]:
    """ Returns for each group of enrichment step, its path to its associated reactions.tsv file.

    Returns
    -------
    Dict[str, str]
        Dictionary associating for each group of enrichment step, its path to its associated reactions.tsv file
    """
    enrich_input_dir = os.path.join(INPUT, ENRICH_D)
    enrich_dict = dict()
    groups_set = get_enrich_groups()
    for group in groups_set:
        enrich_dict[group] = os.path.join(enrich_input_dir, group, REACTIONS_TSV)
    return enrich_dict

def get_enrich_rxn():
    enrich_rxn = set()
    enrich_output_dir = os.path.join(OUTPUT, ENRICH_D)
    for group in os.listdir(enrich_output_dir):
        group_dir = os.path.join(enrich_output_dir, group)
        for file in os.listdir(group_dir):
            if file.endswith('res_validation_networks.tsv'):
                res_file = os.path.join(group_dir, file)
                with open(res_file, 'r') as f:
                    f.__next__()
                    for l in f:
                        enrich_rxn.add(l.split('\t')[0])
    return enrich_rxn


def check_enrich_networks_files(path, ext):
    nw_dir_assoc = {SBML_EXT: SBML_D,
                    PADMET_EXT: PADMET_D}
    nw_dir = os.path.join(path, nw_dir_assoc[ext])
    if not os.path.exists(nw_dir):
        files = [x for x in os.listdir(path) if x.endswith(ext)]
        if len(files) > 0:
            os.mkdir(nw_dir)
            for nw in files:
                os.rename(os.path.join(path, nw), os.path.join(nw_dir, nw))
            return True
        return False
    else:
        files = [x for x in os.listdir(nw_dir) if x.endswith(ext)]
        return len(files) > 0


# CHECK ================================================================================================================

def check_required_files():
    """ Checks if all the required files for the workflow are presents in the good folder

    Raises
    ------
    OSError
        If one of the file is not found
    """
    files_required = [os.path.join(INPUT, TARGETS_D, IN_TARGETS),
                      os.path.join(INPUT, SEEDS_D, IN_SEEDS),
                      os.path.join(INPUT, SEEDS_D, IN_ARTEFACTS),
                      get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
                      ]

    logging.info('Running Check step\n'
                 '==================\n')

    for file in files_required:
        if not os.path.exists(file):
            raise OSError(f'No file {file} found')

    padmet_network = get_file_from_ext(os.path.join(INPUT, NETWORK_D), PADMET_EXT)
    if padmet_network is None:
        sbml_network = get_file_from_ext(os.path.join(INPUT, NETWORK_D), SBML_EXT)
        if sbml_network is None:
            raise OSError(f'No PADMET or SBML network found in {os.path.join(INPUT, NETWORK_D)}')

    check_step_required_files(BLASTP)
    check_step_required_files(ENRICH)

    logging.info('All files required found\n\n---------------\nCheck step done\n')


def check_step_required_files(step: str, group=None) -> bool:
    """ Checks if all the files to run a step are found.

    Parameters
    ----------
    step: str
        Step
    group: str (optional, default=None)
        Specified group for enrichment step

    Returns
    -------
    bool
        True if all required files to run the step are found, False otherwise
    """
    files_step = {BLASTP: [get_file_from_ext(os.path.join(INPUT, DATABASE_D), FASTA_EXT),
                           get_file_from_ext(os.path.join(INPUT, SPECIES_D), FAA_EXT)],
                  ENRICH: get_enrich_reactions_files()}

    if step == BLASTP:
        for file in files_step[step]:
            if file is None:
                logging.info(f'Not all files found to run the {step} step, passing the step\n')
                return False
            elif not os.path.exists(file):
                logging.info(f'Not all files found to run the {step} step, passing the step\n')
                return False
        return True

    elif step == ENRICH:
        if files_step[step] == dict():
            logging.info(f'No group directories for {step} step, passing the step\n')
            return False
        elif group is None:
            all_pres = True
            for g, path in files_step[step].items():
                if not os.path.exists(path):
                    logging.info(f'No reaction file for group {g}, checking PADMET network presence')
                    if not check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, g), PADMET_EXT):
                        logging.info(f'No PADMET network(s) for group {g}, checking SBML network presence')
                        if not check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, g), SBML_EXT):
                            logging.info(f'No reaction file or PADMET networks or SBML networks for group {g}, '
                                         f'--enrich={g} impossible')
                            all_pres = False
                        else:
                            logging.info(f'SBML network(s) found for group {g}, --enrich={g} possible')
                    else:
                        logging.info(f'PADMET network(s) found for group {g}, --enrich={g} possible')
                else:
                    logging.info(f'Reaction file found for group {g}, --enrich={g} possible')
            return all_pres
        else:
            if group != GROUP_ALL:
                if not os.path.exists(files_step[step][group]):
                    logging.info(f'No reaction file for group {group}, checking PADMET network presence')
                    if not check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, group), PADMET_EXT):
                        logging.info(f'No PADMET network for group {group}, checking SBML network presence')
                        if not check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, group), SBML_EXT):
                            logging.info(f'No reaction file or PADMET networks or SBML networks for group {group}, '
                                         f'--enrich={group} impossible')
                            return False
        return True
