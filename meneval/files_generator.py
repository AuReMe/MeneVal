from meneval.environment import *
import os
import logging
from padmet.utils.connection.sbmlGenerator import compound_to_sbml, padmet_to_sbml
from padmet.utils.connection.sbml_to_padmet import sbml_to_padmetSpec
from padmet.utils.exploration.compare_padmet import compare_padmet


def check_file_creation(file: str):
    """ Checks if the file is created (if the file path exists) if not raises error

    Parameters
    ----------
    file: str
        Path of the file to check

    Raises
    -------
    FileNotFoundError
        If the path of the file is not found
    """
    if os.path.exists(file):
        logging.info(f'{file} created')
    else:
        logging.info(f'/!\\ {file} not created !')
        raise FileNotFoundError(f'/!\\ {file} not created !')


def files_exist(files: str or List[str]) -> bool:
    """ Check if a file or a list of files exist.

    Parameters
    ----------
    files: str or list[str]
        File or list of files to check
    Returns
    -------
    bool
        True if all files exist, False if at least one file doesn't exist.
    """
    if type(files) == str:
        files = [files]
    for file in files:
        if not os.path.exists(file):
            return False
    logging.info(f'Files : {files} already exist, passing creation.')
    return True


def generate_targets():
    """ In Input/Targets directory, from the file targets.tsv creates the 3 files :
        - biomass.tsv
        - temp_targets.tsv
        _ targets.sbml
    """
    compartment = 'c'
    if not files_exist([TARGETS[TSV_D], BIOMASS_TSV]):
        with open(os.path.join(INPUT, TARGETS_D, IN_TARGETS), 'r') as in_file, \
                open(TARGETS[TSV_D], 'w') as temp_f, \
                open(BIOMASS_TSV, 'w') as biomass_f:

            # Write biomass reaction information in the file biomass.tsv
            biomass_f.write('reaction_id\tbiomass\n'
                            'comment\treaction of biomass\n'
                            'reversible\tfalse\n'
                            'linked_gene\t\n'
                            '#reactant/product\t#stoichio:compound_id:compar\n')

            # Fill biomass.tsv and temp_targets.tsv with metabolites information
            for line in in_file:
                line = line.strip().split('\t')
                stoichiometry = line[0]
                metabolite = line[1]
                biomass_f.write(f'reactant\t{stoichiometry}:{metabolite}:{compartment}\n')
                temp_f.write(f'{metabolite}\t{compartment}\n')

            # Write biomass product in the file biomass.tsv
            biomass_f.write('product\t1.0:Bio:c\n')
            # Write biomass export reaction in the file biomass.tsv
            biomass_f.write('\nreaction_id\tExport_Bio\n'
                            'comment\tAvoid storage\n'
                            'reversible\tfalse\n'
                            'linked_gene\t\n'
                            'reactant\t1.0:Bio:c\n'
                            'product\t1.0:Bio:C-BOUNDARY\n')

    check_file_creation(TARGETS[TSV_D])
    check_file_creation(BIOMASS_TSV)

    # Create the targets.sbml file from temp_targets.tsv file
    if not files_exist(TARGETS[SBML_D]):
        compound_to_sbml(species_compart=TARGETS[TSV_D], output=TARGETS[SBML_D])
    check_file_creation(TARGETS[SBML_D])


def generate_seeds():
    """ In Input/Seeds directory, from seeds.tsv and artefacts.tsv file, creates the 3 files :
        - seeds_medium.tsv
        - seeds_artefacts.tsv
        - seeds_artefacts.sbml
    """

    seeds_compartment = 'C-BOUNDARY'
    artefacts_compartment = 'c'
    if not files_exist([SEEDS_TSV, SEEDS_ARTEFACTS[TSV_D]]):
        with open(os.path.join(INPUT, SEEDS_D, IN_SEEDS), 'r') as in_seeds, \
                open(os.path.join(INPUT, SEEDS_D, IN_ARTEFACTS), 'r') as in_art, \
                open(SEEDS_TSV, 'w') as seeds_med, \
                open(SEEDS_ARTEFACTS[TSV_D], 'w') as seeds_art:

            # Write seeds metabolites and compartments in seeds_medium.tsv and seeds_artefacts.tsv files
            for line in in_seeds:
                metabolite = line.strip()
                seeds_med.write(f'{metabolite}\t{seeds_compartment}\n')
                seeds_art.write(f'{metabolite}\t{seeds_compartment}\n')

            # Write artefacts metabolites and compartments in seeds_artefacts.tsv file
            for line in in_art:
                metabolite = line.strip()
                seeds_art.write(f'{metabolite}\t{artefacts_compartment}\n')

    check_file_creation(SEEDS_TSV)
    check_file_creation(SEEDS_ARTEFACTS[TSV_D])

    # Create the seeds_artefacts.sbml file from seeds_artefacts.tsv file
    if not files_exist(SEEDS_ARTEFACTS[SBML_D]):
        compound_to_sbml(species_compart=SEEDS_ARTEFACTS[TSV_D], output=SEEDS_ARTEFACTS[SBML_D])
    check_file_creation(SEEDS_ARTEFACTS[SBML_D])


def generate_db_sbml():
    """ In Input/DataBase directory from .padmet dataBase file, create file :
            - database.sbml
    """
    if not files_exist(DB_SBML):
        padmet_to_sbml(padmet=get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT), output=DB_SBML)
    check_file_creation(DB_SBML)


def generate_base_networks():
    """ In Output/Networks directory, create files :
            - 0_base.padmet
            - 0_base.sbml
        """
    # Create 0_medium.padmet
    input_padmet = get_file_from_ext(os.path.join(INPUT, NETWORK_D), PADMET_EXT)
    if input_padmet is None:
        input_sbml = get_file_from_ext(os.path.join(INPUT, NETWORK_D), SBML_EXT)
        input_padmet = input_sbml.split('.')[0] + PADMET_EXT
        sbml_to_padmetSpec(sbml=input_sbml, padmetSpec_file=input_padmet)
        check_file_creation(input_sbml)
    if not files_exist(MEDIUM_NW):
        os.system(f'padmet padmet_medium --padmetSpec={input_padmet} '
                  f'--seeds={SEEDS_TSV}  '
                  f'--padmetRef={get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)} '
                  f'--output={MEDIUM_NW}')
    check_file_creation(MEDIUM_NW)

    # Create 0_base.padmet
    if not files_exist(BASE_NW[PADMET_D]):
        os.system(f'padmet manual_curation --padmetSpec={MEDIUM_NW} '
                  f'--data={BIOMASS_TSV} '
                  f'--output={BASE_NW[PADMET_D]} '
                  f'--category=MANUAL')
    check_file_creation(BASE_NW[PADMET_D])
    # Delete 0_medium.padmet
    os.remove(MEDIUM_NW)

    # Create 0_base.sbml
    if not files_exist(BASE_NW[SBML_D]):
        padmet_to_sbml(padmet=BASE_NW[PADMET_D], output=BASE_NW[SBML_D])
    check_file_creation(BASE_NW[SBML_D])


def generate_reactions_files():
    rxn_files = get_enrich_reactions_files()
    for group, r_path in rxn_files.items():
        if not files_exist(r_path):
            group_path = os.path.join(INPUT, ENRICH_D, group)

            if check_enrich_networks_files(group_path, PADMET_EXT):
                padmet_dir = os.path.join(group_path, PADMET_D)
                compare_padmet(padmet_path=padmet_dir, output=group_path)
                os.remove(os.path.join(group_path, 'metabolites.tsv'))
                os.remove(os.path.join(group_path, 'genes.tsv'))
                os.remove(os.path.join(group_path, 'pathways.tsv'))
                check_file_creation(r_path)

            elif check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, group), SBML_EXT):
                sbml_dir = os.path.join(group_path, SBML_D)
                padmet_dir = os.path.join(group_path, PADMET_D)
                if not os.path.exists(padmet_dir):
                    os.mkdir(padmet_dir)
                for sbml_file in os.listdir(sbml_dir):
                    output_padmet = os.path.join(padmet_dir, str(sbml_file).split('.')[0] + PADMET_EXT)
                    input_sbml = os.path.join(sbml_dir, sbml_file)
                    sbml_to_padmetSpec(sbml=input_sbml, padmetSpec_file=output_padmet)
                check_enrich_networks_files(os.path.join(INPUT, ENRICH_D, group), PADMET_EXT)

                compare_padmet(padmet_path=padmet_dir, output=group_path)
                os.remove(os.path.join(group_path, 'metabolites.tsv'))
                os.remove(os.path.join(group_path, 'genes.tsv'))
                os.remove(os.path.join(group_path, 'pathways.tsv'))
                check_file_creation(r_path)
        else:
            check_file_creation(r_path)
