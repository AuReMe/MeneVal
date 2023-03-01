# -*- coding: utf-8 -*-
"""
Description: Filter the reactions found in solution by Meneco to only keep reactions that lead to a Blast match of the
    proteins associated with the reaction in MetaCyc database against the proteome (or the genome in option).

    Pipeline steps :
    1 - Extract the list of the union of reactions found in solution by Meneco from the file meneco_output.tsv.
    2 - Create a dictionary associating to each reaction the set of Uniprot IDs associated to it in the MetaCyc padmet
        xrefs.
    3 - For each Uniprot ID :
        3.1 - Extract its associated protein sequence from the protein-seq-ids-reduced-70.fasta file.
        3.2 - Perform a Blastp of this sequence against the species proteome and if there is a match write the result to
            the blast_results.tsv file.
        3.3 - If there was no match and the Tblastn option is True: execute a tblastn of this sequence against the
            genome of the species and if there is a match write the result in the blast_results.tsv file.
        3.4 - If there was a match during the Blastp or Tblastn: add the reaction associated with the Uniprot ID to the
            set of reactions to be kept.
    4. Create a new meneco_output_filtered.tsv file retaining only the reaction set that led to a blast match.

::    usage:
        python meneco_validation_BlastP.py -m MENECO_TSV -o OUTPUT_DIR -db METACYC -dbf PROT_FASTA -p PROTEOME [-g GENOME]
    options:
        -m : meneco_output.tsv file : file obtained after execution of padmet enhance_meneco_output command from meneco
            output
        -o : Output directory
        -db : MetaCyc padmet database. Padmet creation from PGDB must have been created with "--prot-ids70" option in
            pgdb_to_padmet --> will assign UNIPROT and PID ids to reactions in their xrefs (with "_70" suffix)
        -dbf : protein-seq-ids-reduced-70.fasta file created while generating MetaCyc padmet file with "--prot-ids70"
            option
        -p : Species proteome
        -g : (optional) Species Genome,  if passed as argument will run TBlastn against species genome if there is no
            match found from the Blastp against species proteome
"""
import os
import logging
import argparse
from padmet.classes.padmetSpec import PadmetSpec
from Bio.Blast.Applications import NcbiblastpCommandline, NcbitblastnCommandline
from Bio import SeqIO
from typing import List, Set


# CONSTANTS ============================================================================================================

# Command line arguments
def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', type=str, required=True, metavar='meneco output file', help='Meneco output file')
    parser.add_argument('-o', type=str, required=True, metavar='directory', help='The output directory')
    parser.add_argument('-db', type=str, required=True, metavar='DataBase padmet', help='DataBase padmet')
    parser.add_argument('-dbf', type=str, required=True, metavar='Proteins Fasta', help='Proteins Fasta')
    parser.add_argument('-p', type=str, required=True, metavar='Proteome', help='Proteome to run BlastP')
    parser.add_argument('-g', type=str, required=False, metavar='Genome', help='Genome to run TBlastn')
    args = parser.parse_args()
    return args.m, args.o, args.db, args.dbf, args.p, args.g


MENECO_TSV, OUTPUT, DB_PADMET, PROT_FASTA, SPECIES_PROTEOME, SPECIES_GENOME = get_command_line_args()


# Parse Files
PADMET_SPEC = PadmetSpec(DB_PADMET)
PROT_FASTA = SeqIO.to_dict(SeqIO.parse(PROT_FASTA, 'fasta'))

TBLASTN = SPECIES_GENOME is not None


# Directories to create
SEQ_DIR = os.path.join(OUTPUT, 'sequences')
RES_DIR = os.path.join(OUTPUT, 'results')

# Blast
BLAST_HEADER = ["sseqid", "evalue", "bitscore", "pident", "length"]
BLAST_RES_FILE = os.path.join(RES_DIR, 'blast_results.tsv')

RXN_PROT_FILE = os.path.join(RES_DIR, 'rxn_prot.tsv')

# Init logger
LOG_FILE = os.path.join(OUTPUT, 'meneco_validation.log')
logging.basicConfig(filename=LOG_FILE, level=logging.INFO, format='%(message)s')
# logging.basicConfig(filename=LOG_FILE, encoding='utf-8', level=logging.INFO, format='%(message)s') python version 3.9


# FUNCTIONS ============================================================================================================

def extract_rxn_from_meneco() -> List[str]:
    """Extract the list of reactions corresponding to the union of reactions from solutions found by Meneco.
    The extraction is done from a meneco_output.tsv file created with padmet enhance_meneco_output.

    Returns
    -------
    List[str]
        List of reactions
    """
    logging.info('Extracting reactions from Meneco output')
    rxn_list = list()
    with open(MENECO_TSV, 'r') as meneco:
        for line in meneco:
            rxn = line.split('\t')[0]
            rxn_list.append(rxn)
    del rxn_list[0]
    logging.info(f'Total of {len(rxn_list)} reactions\n')
    return rxn_list


def get_uniprot_ids_from_rxn(rxn: str) -> Set[str]:
    """ Obtain the set of Uniprot IDs linked to a reaction.

    Parameters
    ----------
    rxn : str
        Reaction to find Uniprot IDs associated with

    Returns
    -------
    Set[str]
        Set of Uniprot IDs liked to the reaction
    """
    logging.info(f'Extract Uniprot IDs linked to reaction {rxn}')
    prot_keys = ('UNIPROT', 'PID')
    suffix = '_70'
    uniprot_set = set()

    try:
        xrefs = PADMET_SPEC.dicOfNode[rxn + '_xrefs'].misc
    except KeyError:
        logging.info(f'0 Protein IDs linked to reaction {rxn}\n')
        return set()

    for prot_key in prot_keys:
        try:
            uniprot_ids = xrefs[prot_key + suffix]
            uniprot_set = uniprot_set.union(set([f'{prot_key}:{prot_id}' for prot_id in uniprot_ids]))
        except KeyError:
            pass

    logging.info(f'{len(uniprot_set)} Protein IDs linked to reaction {rxn}\n')
    return uniprot_set


def get_uniprot_seq(uni_id: str) -> str:
    """Obtain the sequence associated with an Uniprot IDs from the line of the protein-seq-ids-reduced-70.fasta file.

    Parameters
    ----------
    uni_id : str
        Uniprot ID to find the protein sequence associated with.

    Returns
    -------
    str
        The sequence associated with the Uniprot ID.
    """
    logging.info(f'\nGet sequence for {uni_id} protein')
    try:
        seq = PROT_FASTA[uni_id].seq
        return seq
    except KeyError:
        logging.warning(f'\nNo sequence corresponding to {uni_id} in {PROT_FASTA} file.')
        return ''


def create_dirs_and_init_result_file():
    """ Create the 'sequences' and 'results' directories if they don't exist and create the results file."""
    # Directories creation
    if not os.path.exists(SEQ_DIR):
        os.mkdir(SEQ_DIR)
    if not os.path.exists(RES_DIR):
        os.mkdir(RES_DIR)

    # Results file creation and initialization
    f = open(BLAST_RES_FILE, 'w')
    f.write('\t'.join(['Reaction', 'Uniprot ID', 'Sequence', 'E value', 'Bit score', 'Identity (%)', 'Length',
                       'Blast method\n']))
    f.close()
    logging.info(f'Blast results tsv file created in : {BLAST_RES_FILE}')


def run_blastp(prot_seq: str, uniprot_id: str, rxn: str) -> bool:
    """Runs a Blastp of the protein sequence of the Uniprot ID associated with a reaction, against rhe species proteome.
    If there is a match, writes the result in the results file and returns True, otherwise just returns False.

    Parameters
    ----------
    prot_seq : str
        Sequence to use as query for Blastp
    uniprot_id : str
        Uniprot ID corresponding to the protein sequence
    rxn : str
        Reaction within the Uniprot ID has been associated to.

    Returns
    -------
    bool
        True if the Blastp of the protein sequence against the species proteome found a match.
    """
    uniprot_id = uniprot_id.split(':')[1]

    logging.info(f'Blastp of protein {uniprot_id} for reaction {rxn} against {SPECIES_PROTEOME} proteome')

    out_fmt_arg = '"%s %s"' % (6, " ".join(BLAST_HEADER))
    query_fasta = os.path.join(SEQ_DIR, f'{uniprot_id}.fasta')
    if not os.path.exists(query_fasta):
        with open(query_fasta, 'w') as fquery:
            fquery.write(f'>{uniprot_id}\n{prot_seq}\n')

    blastp_output = NcbiblastpCommandline(query=query_fasta, subject=SPECIES_PROTEOME, evalue=1e-10,
                                          outfmt=out_fmt_arg)()[0]
    if blastp_output != '':
        with open(BLAST_RES_FILE, 'a') as res_file:
            output = blastp_output.split('\n')
            for match in output[:-1]:
                res_file.write('\t'.join([rxn, uniprot_id, match, 'Blastp\n']))
        return True
    return False


def run_tblastn(prot_seq: str, uniprot_id: str, rxn: str) -> bool:
    """Runs a Tblastn of the protein sequence of the Uniprot ID associated with a reaction, against rhe species genome.
    If there is a match, writes the result in the results file and returns True, otherwise just returns False.

    Parameters
    ----------
    prot_seq : str
        Sequence to use as query for Tblastn
    uniprot_id : str
        Uniprot ID corresponding to the protein sequence
    rxn : str
        Reaction within the Uniprot ID has been associated to.

    Returns
    -------
    bool
        True if the Tblastn of the protein sequence against the species genome found a match.
    """
    uniprot_id = uniprot_id.split(':')[1]
    logging.info(f'TBlastn of protein {uniprot_id} for reaction {rxn} against {SPECIES_GENOME} genome.')
    out_fmt_arg = '"%s %s"' % (6, " ".join(BLAST_HEADER))
    query_fasta = os.path.join(SEQ_DIR, f'{uniprot_id}.fasta')
    if not os.path.exists(query_fasta):
        with open(query_fasta, 'w') as fquery:
            fquery.write(f'>{uniprot_id}\n{prot_seq}\n')

    tblastn_output = NcbitblastnCommandline(query=query_fasta, subject=SPECIES_GENOME, evalue=1e-10,
                                            outfmt=out_fmt_arg)()[0]
    if tblastn_output != '':
        with open(BLAST_RES_FILE, 'a') as res_file:
            output = tblastn_output.split('\n')
            for match in output[:-1]:
                res_file.write('\t'.join([rxn, uniprot_id, match, 'TBlastn\n']))
        return True
    return False


def create_new_meneco_tsv(kept_rxn: Set[str]):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    new_meneco = os.path.join(RES_DIR, 'meneco_output_filtered.tsv')
    with open(MENECO_TSV, 'r') as in_meneco, open(new_meneco, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                line[-2] = 'Gapfilling BlastP hit'
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {new_meneco}')


def create_rxn_prot_tsv(rxn_prot_dict):
    with open(RXN_PROT_FILE, 'w') as f:
        f.write('Rxn ID\tNb prot IDs\tProt IDs (sep=;)\n')
        for rxn, pid_set in rxn_prot_dict.items():
            f.write(f'{rxn}\t{len(pid_set)}\t{";".join(pid_set)}\n')


# MAIN FUNCTION ========================================================================================================

def meneco_validation():
    """Runs the Meneco validation pipeline :
    1 - Extract the list of the union of reactions found in solution by Meneco from the file meneco_output.tsv.
    2 - Create a dictionary associating to each reaction the set of Uniprot IDs associated to it in the
        protein-seq-ids-reduced-70.dat file.
    3 - For each Uniprot ID :
        3.1 - Extract its associated protein sequence from the protein-seq-ids-reduced-70.seq file.
        3.2 - Perform a Blastp of this sequence against the species proteome and if there is a match write the result to
            the blast_results.tsv file.
        3.3 - If there was no match and the Tblastn option is True: execute a tblastn of this sequence against the
            genome of the species and if there is a match write the result in the blast_results.tsv file.
        3.4 - If there was a match during the Blastp or Tblastn: add the reaction associated with the Uniprot ID to the
            set of reactions to be kept.
    4. Create a new meneco_output_filtered.tsv file retaining only the reaction set that led to a blast match.
    """
    logging.info('Start searching alignments from Meneco output\n'
                 '=============================================\n')

    # 1 - Obtain list of union of reactions retrieved by meneco
    rxn_list = extract_rxn_from_meneco()

    # 2 - Get uniprot_ids linked to all these reactions
    logging.info('Start getting proteins linked to each reaction\n'
                 '==============================================\n')

    uniprot_dict = dict()
    for rxn in rxn_list:
        uniprot_dict[rxn] = get_uniprot_ids_from_rxn(rxn)
    alignment_total = sum(len(uni_id) for uni_id in uniprot_dict.values())
    create_dirs_and_init_result_file()
    create_rxn_prot_tsv(uniprot_dict)

    # Loop on rxn and uniprot IDs
    logging.info('Start blast alignments\n'
                 '======================\n')

    kept_rxn_set = set()
    alignment_number = 0
    for rxn, uniprot_ids in uniprot_dict.items():
        for uni_id in uniprot_ids:
            # Get sequence
            seq = get_uniprot_seq(uni_id)

            # Blastp sp proteome
            keep_rxn = run_blastp(seq, uni_id, rxn)
            # Keep reactions with alignment result
            if keep_rxn:
                kept_rxn_set.add(rxn)

            # Tblastn if not kept
            elif TBLASTN:
                keep_rxn = run_tblastn(seq, uni_id, rxn)
                # Keep reactions with alignment result
                if keep_rxn:
                    kept_rxn_set.add(rxn)

            # Inform advancement
            alignment_number += 1
            logging.info(f'{alignment_number}/{alignment_total} done')

    logging.info(f'\n\n{len(kept_rxn_set)} reactions with alignment : {", ".join(kept_rxn_set)}')
    create_new_meneco_tsv(kept_rxn_set)
    logging.info(f'=================\n'
                 f'End of alignments')


# MAIN SCRIPT ==========================================================================================================

meneco_validation()
