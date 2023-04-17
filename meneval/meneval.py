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
from files_generator import *
from meneco_utils import *
from stats_recap import *
from validation_BlastP import validation_blastp
from validation_networks import validation_networks
import os
import shutil
import argparse
import logging


# ARGUMENTS
def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init', action='store_true', required=False, help='Init environment')
    parser.add_argument('--check', action='store_true', required=False, help='Check for required input files')
    parser.add_argument('--files', action='store_true', required=False, help='generate additional required input files')
    parser.add_argument('--blastp', action='store_true', required=False, help='Runs blastp step')
    parser.add_argument('--holobiont', action='store_true', required=False, help='Runs holobiont step')
    parser.add_argument('--aucome', action='store_true', required=False, help='Runs Aucome step')
    parser.add_argument('--group', type=str, required=False, metavar='group name', help='Group name for aucome step')
    parser.add_argument('--fill', action='store_true', required=False, help='Runs fill step')
    parser.add_argument('--workflow', action='store_true', required=False, help='Runs all steps')
    args = parser.parse_args()
    return args.init, args.check, args.files, args.blastp, args.holobiont, args.aucome, args.group, args.fill, \
        args.workflow


INIT, CHECK, FILES_GENERATION, BLASTP_STEP, HOLOBIONT_STEP, AUCOME_STEP, GROUP, FILL_STEP, WORKFLOW = \
    get_command_line_args()

if WORKFLOW:
    INIT = False
    CHECK = True
    FILES_GENERATION = True
    BLASTP_STEP = True
    HOLOBIONT_STEP = True
    AUCOME_STEP = True
    FILL_STEP = True


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
    # Run function
    validation_networks(name, output, meneco_tsv, reactions_tsv)
    # Move filtered tsv output file in correct path
    os.rename(os.path.join(output, 'meneco_output_filtered.tsv'),
              meneco_filtered)


# 3RD VALIDATION
def aucome_step(meneco_tsv, meneco_filtered):
    name = GROUP
    output = os.path.join(OUTPUT, AUCOME_D)
    reactions_tsv = os.path.join(INPUT, AUCOME_D, REACTIONS_TSV)
    group_template = os.path.join(INPUT, AUCOME_D, GROUPS_TSV)
    # Run function
    validation_networks(name, output, meneco_tsv, reactions_tsv, group_template, GROUP)
    # Move filtered tsv output file in correct path
    os.rename(os.path.join(output, 'meneco_output_filtered.tsv'),
              meneco_filtered)


def final_step(meneco_tsv, meneco_filtered):
    shutil.copy(meneco_tsv, meneco_filtered)


def run_step(num):
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
            raise OSError(f'No padmet and sbml files found for running step {name}')

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

    # TSV output creation
    logging.info(f'\nCreate Meneco tsv output :\n{40 * "-"}\n')
    if not os.path.exists(meneco_tsv):
        meneco_out_txt_to_tsv(meneco_out, meneco_tsv)
    else:
        logging.info(f'{meneco_tsv} file found, passing tsv creation.')

    # Validation step
    logging.info(f'\nRunning {name} validation step :\n{40 * "-"}\n')
    if not os.path.exists(meneco_filtered):
        if name == names_list[0]:
            blastp_step(meneco_tsv, meneco_filtered)
            add_genes_tsv(meneco_filtered)
        if name == names_list[1]:
            holobiont_step(meneco_tsv, meneco_filtered)
        if name == names_list[2]:
            aucome_step(meneco_tsv, meneco_filtered)
        if name == names_list[3]:
            final_step(meneco_tsv, meneco_filtered)
    else:
        logging.info(f'{meneco_filtered} file found, passing {name} validation.')

    # Add reactions to network
    logging.info(f'\nAdding reactions found to network :\n{40 * "-"}\n')
    if not os.path.exists(dict_nw[PADMET_D]):
        add_rxn_to_nw(prev_network_padmet, dict_nw[PADMET_D], meneco_filtered)
    else:
        logging.info(f'{dict_nw[PADMET_D]} file found, passing adding reactions to network.')

    # Create SBML network
    logging.info(f'\nConvert Padmet to SBML:\n{40 * "-"}\n')
    if not os.path.exists(dict_nw[SBML_D]):
        padmet_to_sbml(padmet=dict_nw[PADMET_D], output=dict_nw[SBML_D])
    else:
        logging.info(f'{dict_nw[SBML_D]} file found, passing convert Padmet to SBML.')

    logging.info(f'\n-----------\nStep {num} Done\n')


# INITIALIZATION AND CHECK =============================================================================================
if INIT:
    create_folders()
if CHECK:
    check_required_files()


# GENERATE MENECO FILES NEEDED =========================================================================================
if FILES_GENERATION:
    generate_seeds()
    generate_targets()
    generate_db_sbml()
    generate_base_networks()
    logging.info('All files needed created successfully')


# RUN STEPS ============================================================================================================

if BLASTP_STEP:
    run_step(1)

if HOLOBIONT_STEP:
    run_step(2)

if AUCOME_STEP:
    run_step(3)

if FILL_STEP:
    run_step(4)
    make_meneco_stats()
