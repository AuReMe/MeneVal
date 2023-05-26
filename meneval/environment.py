"""
Functions to initialize and check the input files and folder required for meneco_validation workflow.
The initialisation will create all required folders.
The check will verify if all the required files are placed in the corrects folders.
"""
import os
import logging
from typing import Dict, List

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

MEDIUM_NW = os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'1_medium{PADMET_EXT}')
BASE_NW = {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'1_base{PADMET_EXT}'),
           SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'1_base{SBML_EXT}')}

# # BLASTP GAPFILLING NETWORKS
# BLASTP_GF_NW = {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'2_gapfilling_blastp{PADMET_EXT}'),
#                 SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'2_gapfilling_blastp{SBML_EXT}')}
#
# # ENRICHMENT GAPFILLING NETWORKS
# HOLOBIONT_GF_NW = {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'3_gapfilling_holobiont{PADMET_EXT}'),
#                    SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'3_gapfilling_holobiont{SBML_EXT}')}
#
#
# # FINAL GAPFILLING NETWORKS
# FINAL_GF_NW = {PADMET_D: os.path.join(OUTPUT, NETWORK_D, PADMET_D, f'5_gapfilling_final{PADMET_EXT}'),
#                SBML_D: os.path.join(OUTPUT, NETWORK_D, SBML_D, f'5_gapfilling_final{SBML_EXT}')}


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


def get_enrich_reactions_files():
    enrich_input_dir = os.path.join(INPUT, ENRICH_D)
    enrich_dict = dict()
    if os.listdir(enrich_input_dir) != list():
        for directory in os.listdir(enrich_input_dir):
            enrich_dict[directory] = os.path.join(enrich_input_dir, directory, REACTIONS_TSV)
    else:
        print(f'Y a rien dans {enrich_input_dir}')
    return enrich_dict


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
                      get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT),
                      get_file_from_ext(os.path.join(INPUT, NETWORK_D), PADMET_EXT)
                      ]

    logging.info('Running Check step\n'
                 '==================\n')

    for file in files_required:
        if not os.path.exists(file):
            raise OSError(f'No file {file} found')

    check_step_required_files(1)
    check_step_required_files(2)
    check_step_required_files(3)

    logging.info('All files required found\n\n---------------\nCheck step done\n')


def check_step_required_files(step_num: int, group=None) -> bool:
    """ Checks if all the files to run a step are found.

    Parameters
    ----------
    step_num: int
        Number of the step (1=blastp, 2=holobiont, 3=aucome)

    Returns
    -------
    bool
        True if all required files to run the step are found, False otherwise
    """
    names_list = ['', 'BLASTP', 'ENRICHMENT']
    files_step = {1: [get_file_from_ext(os.path.join(INPUT, DATABASE_D), FASTA_EXT),
                      get_file_from_ext(os.path.join(INPUT, SPECIES_D), FAA_EXT)],
                  2: get_enrich_reactions_files()}

    if step_num == 1:
        for file in files_step[step_num]:
            if file is None:
                logging.info(f'Not all files found to run the {names_list[step_num]} step, passing the step\n')
                return False
            elif not os.path.exists(file):
                logging.info(f'Not all files found to run the {names_list[step_num]} step, passing the step\n')
                return False
        return True

    elif step_num == 2:
        if files_step[step_num] == dict():
            logging.info(f'Not all files found to run the {names_list[step_num]} step, passing the step\n')
            return False
        elif group is None:
            for group, path in files_step[step_num].items():
                pass
        else:
            pass

        return True



print(check_step_required_files(2))
