"""
Functions to initialize and check the input files and folder required for meneco_validation workflow.
The initialisation will create all required folders.
The check will verify if all the required files are placed in the corrects folders.
"""
import os
from typing import Dict, List

# BASE ENVIRONMENT =====================================================================================================

# MAIN DIRECTORIES
INPUT = os.path.join('Input')
OUTPUT = os.path.join('Output')

AUCOME_D = 'AuCoMe'
HOLOBIONT_D = 'Holobiont'
DATABASE_D = 'DataBase'
NETWORK_D = 'Network'
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
TOOL_OUTPUTS_D = 'Tool_outputs'

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
        print(f'No file with the extension {ext} in path {path}')
    elif len(files) > 1:
        print(f'More than 1 file with the extension {ext} in path {path}')
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
MEDIUM_NW = os.path.join(INPUT, NETWORK_D, PADMET_D, f'1_medium{PADMET_EXT}')
BASE_NW = {PADMET_D: os.path.join(INPUT, NETWORK_D, PADMET_D, f'1_base{PADMET_EXT}'),
           SBML_D: os.path.join(INPUT, NETWORK_D, SBML_D, f'1_base{SBML_EXT}')}

# BLASTP GAPFILLING NETWORKS
BLASTP_GF_NW = {PADMET_D: os.path.join(INPUT, NETWORK_D, PADMET_D, f'2_gapfilling_blastp{PADMET_EXT}'),
                SBML_D: os.path.join(INPUT, NETWORK_D, SBML_D, f'2_gapfilling_blastp{SBML_EXT}')}

# HOLOBIONT GAPFILLING NETWORKS
HOLOBIONT_GF_NW = {PADMET_D: os.path.join(INPUT, NETWORK_D, PADMET_D, f'3_gapfilling_holobiont{PADMET_EXT}'),
                   SBML_D: os.path.join(INPUT, NETWORK_D, SBML_D, f'3_gapfilling_holobiont{SBML_EXT}')}

# AUCOME GAPFILLING NETWORKS
AUCOME_GF_NW = {PADMET_D: os.path.join(INPUT, NETWORK_D, PADMET_D, f'4_gapfilling_aucome{PADMET_EXT}'),
                SBML_D: os.path.join(INPUT, NETWORK_D, SBML_D, f'4_gapfilling_aucome{SBML_EXT}')}

# FINAL GAPFILLING NETWORKS
FINAL_GF_NW = {PADMET_D: os.path.join(INPUT, NETWORK_D, PADMET_D, f'5_gapfilling_final{PADMET_EXT}'),
               SBML_D: os.path.join(INPUT, NETWORK_D, SBML_D, f'5_gapfilling_final{SBML_EXT}')}


# INITIALIZATION =======================================================================================================

def create_folders():
    """ Creates all the input and outputs folders required for the workflow
    """
    dir_archi = {INPUT: [AUCOME_D,
                         DATABASE_D,
                         HOLOBIONT_D,
                         SPECIES_D,
                         {NETWORK_D: [PADMET_D,
                                      SBML_D]},
                         SEEDS_D,
                         TARGETS_D,
                         ],

                 OUTPUT: [AUCOME_D,
                          BLASTP_D,
                          HOLOBIONT_D,
                          {MENECO_D: [FILTERED_D,
                                      TOOL_OUTPUTS_D,
                                      TSV_D]}
                          ]
                 }
    create_dir_rec(dir_archi)


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
            print(f'{parent_path} created')
        else:
            print(f'{parent_path} already exists')
        for child in child_list:
            if type(child) == str:
                child_path = os.path.join(path, parent_rep, child)
                if not os.path.exists(child_path):
                    os.mkdir(child_path)
                    print(f'{child_path} created')
                else:
                    print(f'{child_path} already exists')
            else:
                create_dir_rec(child, os.path.join(path, parent_rep))


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

    files_blastp = [get_file_from_ext(os.path.join(INPUT, DATABASE_D), FASTA_EXT),
                    get_file_from_ext(os.path.join(INPUT, SPECIES_D), FAA_EXT)]

    files_holobiont = [os.path.join(INPUT, HOLOBIONT_D, REACTIONS_TSV)]

    files_aucome = [os.path.join(INPUT, AUCOME_D, GROUPS_TSV),
                    os.path.join(INPUT, AUCOME_D, REACTIONS_TSV)]

    for file in files_required:
        if not os.path.exists(file):
            raise OSError(f'No file {file} found')

    for file in files_blastp:
        if not os.path.exists(file):
            print('Not all files found to run the blastp step, add it if you want to run this step')

    for file in files_holobiont:
        if not os.path.exists(file):
            print('Not all files found to run the holobiont step, add it if you want to run this step')

    for file in files_aucome:
        if not os.path.exists(file):
            print('Not all files found to run the aucome step, add it if you want to run this step')

    print('All files required found -> Check passed')
