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
"""
import os
import logging
import padmet.classes.padmetSpec
from padmet.classes.padmetSpec import PadmetSpec
from Bio.Blast.Applications import NcbiblastpCommandline, NcbitblastnCommandline
from Bio import SeqIO
from typing import List, Set


# CONSTANTS ============================================================================================================

# Blast
BLAST_HEADER = ["sseqid", "evalue", "bitscore", "pident", "length"]


# ENVIRONMENT ==========================================================================================================
def get_directories(output):
    seq_dir = os.path.join(output, 'sequences')
    res_dir = os.path.join(output, 'results')

    blast_res_file = os.path.join(res_dir, 'blast_results.tsv')
    rxn_prot_file = os.path.join(res_dir, 'rxn_prot.tsv')

    log_file = os.path.join(output, 'meneco_validation.log')

    return seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file


def init_logger(log_file):
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(message)s')
    # logging.basicConfig(filename=log_file, encoding='utf-8', level=logging.INFO, format='%(message)s') python v 3.9

# FUNCTIONS ============================================================================================================


def extract_rxn_from_meneco(meneco_tsv: str) -> List[str]:
    """Extract the list of reactions corresponding to the union of reactions from solutions found by Meneco.
    The extraction is done from a meneco_output.tsv file created with padmet enhance_meneco_output.

    Returns
    -------
    List[str]
        List of reactions
    """
    logging.info('Extracting reactions from Meneco output')
    rxn_list = list()
    with open(meneco_tsv, 'r') as meneco:
        for line in meneco:
            rxn = line.split('\t')[0]
            rxn_list.append(rxn)
    del rxn_list[0]
    logging.info(f'Total of {len(rxn_list)} reactions\n')
    return rxn_list


def get_uniprot_ids_from_rxn(rxn: str, padmet_spec: padmet.classes.padmetSpec.PadmetSpec) -> Set[str]:
    """ Obtain the set of Uniprot IDs linked to a reaction.

    Parameters
    ----------
    rxn : str
        Reaction to find Uniprot IDs associated with
    padmet_spec : padmet.classes.padmetSpec.PadmetSpec
        PadmetSpec of DataBase

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
        xrefs = padmet_spec.dicOfNode[rxn + '_xrefs'].misc
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


def get_uniprot_seq(uni_id: str, prot_fasta: dict) -> str:
    """Obtain the sequence associated with an Uniprot IDs from the line of the protein-seq-ids-reduced-70.fasta file.

    Parameters
    ----------
    uni_id : str
        Uniprot ID to find the protein sequence associated with.
    prot_fasta: dict
        Dictionary of records of proteins sequences

    Returns
    -------
    str
        The sequence associated with the Uniprot ID.
    """
    logging.info(f'\nGet sequence for {uni_id} protein')
    try:
        seq = prot_fasta[uni_id].seq
        return seq
    except KeyError:
        logging.warning(f'\nNo sequence corresponding to {uni_id} in {prot_fasta} file.')
        return ''


def create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file):
    """ Create the 'sequences' and 'results' directories if they don't exist and create the results file."""
    # Directories creation
    if not os.path.exists(seq_dir):
        os.mkdir(seq_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)

    # Results file creation and initialization
    f = open(blast_res_file, 'w')
    f.write('\t'.join(['Reaction', 'Uniprot ID', 'Sequence', 'E value', 'Bit score', 'Identity (%)', 'Length',
                       'Blast method\n']))
    f.close()
    logging.info(f'Blast results tsv file created in : {blast_res_file}')


def run_blastp(prot_seq: str, uniprot_id: str, rxn: str, species_proteome: str, seq_dir: str, blast_res_file: str) \
        -> bool:
    """Runs a Blastp of the protein sequence of the Uniprot ID associated with a reaction, against rhe species proteome.
    If there is a match, writes the result in the results file and returns True, otherwise just returns False.

    Parameters
    ----------
    prot_seq: str
        Sequence to use as query for Blastp
    uniprot_id: str
        Uniprot ID corresponding to the protein sequence
    rxn: str
        Reaction within the Uniprot ID has been associated to.
    species_proteome: str
        Proteome of the species
    seq_dir: str
        Directory where sequences are stored
    blast_res_file: str
        File to write results information of the blast


    Returns
    -------
    bool
        True if the Blastp of the protein sequence against the species proteome found a match.
    """
    uniprot_id = uniprot_id.split(':')[1]

    logging.info(f'Blastp of protein {uniprot_id} for reaction {rxn} against {species_proteome} proteome')

    out_fmt_arg = '"%s %s"' % (6, " ".join(BLAST_HEADER))
    query_fasta = os.path.join(seq_dir, f'{uniprot_id}.fasta')
    if not os.path.exists(query_fasta):
        with open(query_fasta, 'w') as fquery:
            fquery.write(f'>{uniprot_id}\n{prot_seq}\n')

    blastp_output = NcbiblastpCommandline(query=query_fasta, subject=species_proteome, evalue=1e-10,
                                          outfmt=out_fmt_arg)()[0]
    if blastp_output != '':
        with open(blast_res_file, 'a') as res_file:
            output = blastp_output.split('\n')
            for match in output[:-1]:
                res_file.write('\t'.join([rxn, uniprot_id, match, 'Blastp\n']))
        return True
    return False


def run_tblastn(prot_seq: str, uniprot_id: str, rxn: str, species_genome: str, seq_dir: str, blast_res_file: str) \
        -> bool:
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
    species_genome: str
        Genome of the species
    seq_dir: str
        Directory where sequences are stored
    blast_res_file: str
        File to write results information of the blast

    Returns
    -------
    bool
        True if the Tblastn of the protein sequence against the species genome found a match.
    """
    uniprot_id = uniprot_id.split(':')[1]
    logging.info(f'TBlastN of protein {uniprot_id} for reaction {rxn} against {species_genome} genome.')
    out_fmt_arg = '"%s %s"' % (6, " ".join(BLAST_HEADER))
    query_fasta = os.path.join(seq_dir, f'{uniprot_id}.fasta')
    if not os.path.exists(query_fasta):
        with open(query_fasta, 'w') as fquery:
            fquery.write(f'>{uniprot_id}\n{prot_seq}\n')

    tblastn_output = NcbitblastnCommandline(query=query_fasta, subject=species_genome, evalue=1e-10,
                                            outfmt=out_fmt_arg)()[0]
    if tblastn_output != '':
        with open(blast_res_file, 'a') as res_file:
            output = tblastn_output.split('\n')
            for match in output[:-1]:
                res_file.write('\t'.join([rxn, uniprot_id, match, 'TBlastN\n']))
        return True
    return False


def create_new_meneco_tsv(meneco_tsv: str, kept_rxn: Set[str], res_dir: str):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    meneco_tsv: str
        Meneco results in TSV format (output from enhanced_meneco_output in padmet)
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    res_dir: str
        Directory to store results
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    new_meneco = os.path.join(res_dir, 'meneco_output_filtered.tsv')
    with open(meneco_tsv, 'r') as in_meneco, open(new_meneco, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                line[-2] = 'Gap-filling BlastP hit'
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {new_meneco}')


def create_rxn_prot_tsv(rxn_prot_dict, rxn_prot_file):
    with open(rxn_prot_file, 'w') as f:
        f.write('Rxn ID\tNb prot IDs\tProt IDs (sep=;)\n')
        for rxn, pid_set in rxn_prot_dict.items():
            f.write(f'{rxn}\t{len(pid_set)}\t{";".join(pid_set)}\n')


# MAIN FUNCTION ========================================================================================================

def meneco_validation_blast(meneco_tsv, output, db_padmet, prot_fasta, species_proteome, species_genome=None):
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
    padmet_spec = PadmetSpec(db_padmet)
    prot_fasta = SeqIO.to_dict(SeqIO.parse(prot_fasta, 'fasta'))
    tblastn = species_genome is not None
    seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(output)
    init_logger(log_file)

    logging.info('Start searching alignments from Meneco output\n'
                 '=============================================\n')

    # 1 - Obtain list of union of reactions retrieved by meneco
    rxn_list = extract_rxn_from_meneco(meneco_tsv)

    # 2 - Get uniprot_ids linked to all these reactions
    logging.info('Start getting proteins linked to each reaction\n'
                 '==============================================\n')

    uniprot_dict = dict()
    for rxn in rxn_list:
        uniprot_dict[rxn] = get_uniprot_ids_from_rxn(rxn, padmet_spec)
    alignment_total = sum(len(uni_id) for uni_id in uniprot_dict.values())
    create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file)
    create_rxn_prot_tsv(uniprot_dict, rxn_prot_file)

    # Loop on rxn and uniprot IDs
    logging.info('Start blast alignments\n'
                 '======================\n')

    kept_rxn_set = set()
    alignment_number = 0
    for rxn, uniprot_ids in uniprot_dict.items():
        for uni_id in uniprot_ids:
            # Get sequence
            seq = get_uniprot_seq(uni_id, prot_fasta)

            # Blastp sp proteome
            keep_rxn = run_blastp(seq, uni_id, rxn, species_proteome, seq_dir, blast_res_file)
            # Keep reactions with alignment result
            if keep_rxn:
                kept_rxn_set.add(rxn)

            # Tblastn if not kept
            elif tblastn:
                keep_rxn = run_tblastn(seq, uni_id, rxn, species_genome, seq_dir, blast_res_file)
                # Keep reactions with alignment result
                if keep_rxn:
                    kept_rxn_set.add(rxn)

            # Inform advancement
            alignment_number += 1
            logging.info(f'{alignment_number}/{alignment_total} done')

    logging.info(f'\n\n{len(kept_rxn_set)} reactions with alignment : {", ".join(kept_rxn_set)}')
    create_new_meneco_tsv(meneco_tsv, kept_rxn_set, res_dir)
    logging.info(f'=================\n'
                 f'End of alignments')
