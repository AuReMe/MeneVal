from environment import *
import os
from padmet.utils.connection.sbmlGenerator import compound_to_sbml, padmet_to_sbml


def generate_targets():
    """ In Targets directory, from the file targets.tsv creates the 3 files :
        - biomass.tsv
        - temp_targets.tsv
        _ targets.sbml
    """
    compartment = 'c'
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

    print(f'{TARGETS[TSV_D]} created')
    print(f'{BIOMASS_TSV} created')

    # Create the targets.sbml file from temp_targets.tsv file
    compound_to_sbml(species_compart=TARGETS[TSV_D], output=TARGETS[SBML_D])
    print(f'{TARGETS[SBML_D]} created')


def generate_seeds():
    """ In Seeds directory, from seeds.tsv and artefacts.tsv file, creates the 3 files :
        - seeds_medium.tsv
        - seeds_artefacts.tsv
        - seeds_artefacts.sbml
    """
    seeds_compartment = 'C-BOUNDARY'
    artefacts_compartment = 'c'
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

    print(f'{SEEDS_TSV} created')
    print(f'{SEEDS_ARTEFACTS[TSV_D]} created')

    # Create the seeds_artefacts.sbml file from seeds_artefacts.tsv file
    compound_to_sbml(species_compart=SEEDS_ARTEFACTS[TSV_D], output=SEEDS_ARTEFACTS[SBML_D])
    print(f'{SEEDS_ARTEFACTS[SBML_D]} created')


def generate_db_sbml():
    padmet_to_sbml(padmet=get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT), output=DB_SBML)
    print(f'{DB_SBML} created')


def generate_base_networks():
    """ In Networks directory, files :
            - 1_<SPECIES>_medium.padmet
            - 1_<SPECIES>_biomass.padmet
            - 1_<SPECIES>_base.sbml
        """
    # Create 1_<SPECIES>_medium.padmet
    base_padmet = get_file_from_ext(os.path.join(INPUT, NETWORK_D), PADMET_EXT)
    os.system(f'padmet padmet_medium --padmetSpec={base_padmet} '
              f'--seeds={SEEDS_TSV}  '
              f'--padmetRef={get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)} '
              f'--output={MEDIUM_NW}')
    print(f'{MEDIUM_NW} created')

    # Create 1_<SPECIES>_biomass.padmet
    os.system(f'padmet manual_curation --padmetSpec={MEDIUM_NW} '
              f'--data={BIOMASS_TSV} '
              f'--output={BASE_NW[PADMET_D]} '
              f'--category=MANUAL')
    print(f'{BASE_NW[PADMET_D]} created')

    # Create 1_<SPECIES>_base.sbml
    padmet_to_sbml(padmet=BASE_NW[PADMET_D], output=BASE_NW[SBML_D])
    print(f'{BASE_NW[SBML_D]} created')
