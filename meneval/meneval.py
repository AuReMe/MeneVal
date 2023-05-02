"""
INPUT FILES TO ADD :
    |- Input
        |- AuCoMe
            |- reactions.tsv : File reactions.tsv obtained from analysis step in AuCoMe run
            |- group_template.tsv (opt) : File of group template for group species selection
        |- DataBase
            |- <DataBase>.padmet : DataBase with ".padmet" extension created from PGDB. Must have been created with
                "--prot-ids70" option in pgdb_to_padmet --> will assign UNIPROT and PID ids to reactions in their xrefs
                (with "_70" suffix)
            |- <prot_seq>.fasta : protein-seq-ids-reduced-70.fasta file created while generating MetaCyc padmet file
                with "--prot-ids70" option
        |- Holobiont
            |- reactions.tsv : File reactions.tsv obtained with compare padmet from holobiont networks
            |- group_template.tsv (opt) : File of group template for group species selection
        |- Networks
            |- <sp>.padmet : Networks in ".padmet" for gap filling
        |- Seeds
            |- <seeds>.tsv : Seeds in ".tsv" (artefacts included)
            |- <artefacts>.tsv : Artefacts in .tsv format
        |- Species_seq
            |- <sp>.faa : Species proteome
            |- <sp>.fna (opt) : Species genome if tblastn
"""
# IMPORTS
from meneval.files_generator import *
from meneval.meneco_utils import *
from meneval.stats_recap import *
from meneval.validation_BlastP import validation_blastp
from meneval.validation_networks import validation_networks
import os
import shutil
import logging


# Init logger
LOG_FILE = os.path.join('meneco_validation.log')
logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(message)s')

# FUNCTIONS ============================================================================================================


# 1ST VALIDATION
def blastp_step(meneco_tsv, meneco_filtered):
    # Get files in dirs
    output = os.path.join(OUTPUT, BLASTP_D)
    db_padmet = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    prot_fasta = get_file_from_ext(os.path.join(INPUT, DATABASE_D), FASTA_EXT)
    species_proteome = get_file_from_ext(os.path.join(INPUT, SPECIES_D), FAA_EXT)
    species_genome = get_file_from_ext(os.path.join(INPUT, SPECIES_D), '.fna')
    # Run functions
    rxn_list = extract_rxn_from_meneco(meneco_tsv)
    kept_rxn_set = validation_blastp(rxn_list, output, db_padmet, prot_fasta, species_proteome, species_genome)
    create_new_meneco_tsv(meneco_tsv, kept_rxn_set, meneco_filtered, 'Gap-filling BlastP hit')


# 2ND VALIDATION
def holobiont_step(meneco_tsv, meneco_filtered):
    name = 'Holobiont'
    output = os.path.join(OUTPUT, HOLOBIONT_D)
    reactions_tsv = os.path.join(INPUT, HOLOBIONT_D, REACTIONS_TSV)
    # Run functions
    rxn_list = extract_rxn_from_meneco(meneco_tsv)
    kept_rxn_set = validation_networks(name, output, rxn_list, reactions_tsv)
    create_new_meneco_tsv(meneco_tsv, kept_rxn_set, meneco_filtered, f'Potential {name} source')


# 3RD VALIDATION
def aucome_step(meneco_tsv, meneco_filtered, group):
    name = group
    output = os.path.join(OUTPUT, AUCOME_D)
    reactions_tsv = os.path.join(INPUT, AUCOME_D, REACTIONS_TSV)
    group_template = os.path.join(INPUT, AUCOME_D, GROUPS_TSV)
    # Run function
    rxn_list = extract_rxn_from_meneco(meneco_tsv)
    kept_rxn_set = validation_networks(name, output, rxn_list, reactions_tsv, group_template, group)
    create_new_meneco_tsv(meneco_tsv, kept_rxn_set, meneco_filtered, f'Potential {name} source')


def final_step(meneco_tsv, meneco_filtered):
    shutil.copy(meneco_tsv, meneco_filtered)


def run_step(num, group=None):
    # Get appropriated file
    names_list = ['BLASTP', 'HOLOBIONT', 'AUCOME', 'FILL']
    name = names_list[num - 1]
    networks_rank = [BASE_NW, BLASTP_GF_NW, HOLOBIONT_GF_NW, AUCOME_GF_NW, FINAL_GF_NW]

    dict_nw = networks_rank[num]
    prev = num - 1
    while (not os.path.exists(networks_rank[prev][SBML_D] or
           not os.path.exists(networks_rank[prev][PADMET_D]))):
        prev -= 1
        if prev == -1:
            raise FileNotFoundError(f'No padmet and sbml files found for running step {name}')

    prev_network_sbml = networks_rank[prev][SBML_D]
    prev_network_padmet = networks_rank[prev][PADMET_D]

    # Begin step
    logging.info(f'{50 * "="}\n\tSTEP {num} : MENECO + {name} VALIDATION\n{50 * "="}\n')
    meneco_out, meneco_tsv, meneco_filtered = get_meneco_files(num)

    # Meneco run
    logging.info(f'Running Meneco :\n{40 * "-"}\n')
    if not os.path.exists(meneco_out):
        run_meneco(prev_network_sbml, meneco_out)
    else:
        logging.info(f'{meneco_out} file found, passing meneco run.')
    check_file_creation(meneco_out)

    if exists_not_producible_targets(meneco_out):

        # TSV output creation
        logging.info(f'\nCreate Meneco tsv output :\n{40 * "-"}\n')
        if not os.path.exists(meneco_tsv):
            meneco_json_to_tsv(meneco_out, meneco_tsv)
        else:
            logging.info(f'{meneco_tsv} file found, passing tsv creation.')
        check_file_creation(meneco_tsv)

        # Validation step
        logging.info(f'\nRunning {name} validation step :\n{40 * "-"}\n')
        if not os.path.exists(meneco_filtered):
            if name == names_list[0]:
                blastp_step(meneco_tsv, meneco_filtered)
                add_genes_tsv(meneco_filtered)
            if name == names_list[1]:
                holobiont_step(meneco_tsv, meneco_filtered)
            if name == names_list[2]:
                aucome_step(meneco_tsv, meneco_filtered, group)
            if name == names_list[3]:
                final_step(meneco_tsv, meneco_filtered)
            logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(message)s')
        else:
            logging.info(f'{meneco_filtered} file found, passing {name} validation.')
        check_file_creation(meneco_filtered)

        # Add reactions to network
        logging.info(f'\nAdding reactions found to network :\n{40 * "-"}\n')
        if not os.path.exists(dict_nw[PADMET_D]):
            add_rxn_to_nw(prev_network_padmet, dict_nw[PADMET_D], meneco_filtered)
        else:
            logging.info(f'{dict_nw[PADMET_D]} file found, passing adding reactions to network.')
        check_file_creation(dict_nw[PADMET_D])

        # Create SBML network
        logging.info(f'\nConvert Padmet to SBML:\n{40 * "-"}\n')
        if not os.path.exists(dict_nw[SBML_D]):
            padmet_to_sbml(padmet=dict_nw[PADMET_D], output=dict_nw[SBML_D])
        else:
            logging.info(f'{dict_nw[SBML_D]} file found, passing convert Padmet to SBML.')
        check_file_creation(dict_nw[SBML_D])

    else:
        logging.info('\n--> No targets left to reach. Finishing step.')

    logging.info(f'\n-----------\nStep {num} Done\n')


def generate_files():
    logging.info('Running files generation step :\n'
                 '===============================\n')
    generate_seeds()
    generate_targets()
    generate_db_sbml()
    generate_base_networks()
    logging.info('\nAll files needed created successfully')
    logging.info('\n--------------------------\nFiles generation step done\n')

