import meneco
import json
import logging
import os
from typing import Tuple, Set
from meneval.environment import *


def extract_rxn_from_meneco(meneco_tsv: str) -> List[str]:
    """Extract the list of reactions corresponding to the union of reactions from solutions found by Meneco.
    The extraction is done from a meneco_output.tsv file created with padmet enhance_meneco_output.

    Parameters
    ----------
    meneco_tsv: str
        Meneco results in TSV format (output from enhanced_meneco_output in padmet)

    Returns
    -------
    List[str]
        List of reactions
    """
    logging.info('Extracting reactions from Meneco output')
    rxn_list = list()
    with open(meneco_tsv, 'r') as m_tsv:
        for line in m_tsv:
            rxn = line.split('\t')[0]
            rxn_list.append(rxn)
    del rxn_list[0]
    logging.info(f'Total of {len(rxn_list)} reactions\n')
    return rxn_list


def create_new_meneco_tsv(meneco_tsv: str, kept_rxn: Set[str], output_meneco_tsv: str, message: str or Dict[str, str]):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    meneco_tsv: str
        Meneco results in TSV format (output from enhanced_meneco_output in padmet)
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    output_meneco_tsv: str
        Path to the new meneco tsv output
    message: str or Dict[str, str]
        Comment for adding the reaction
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    with open(meneco_tsv, 'r') as in_meneco, open(output_meneco_tsv, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                if type(message) == str:
                    line[-2] = message
                else:
                    line[-2] = message[line[0]]
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {output_meneco_tsv}')


def get_meneco_files(num: int) -> Tuple[str, str, str]:
    """ Returns the files names of meneco outputs (meneco json output, meneco tsv output and filtered meneco tsv output)
    according the number of the step results desired.

    Parameters
    ----------
    num: int
        Number defining the meneco files to get. (corresponding to a step)

    Returns
    -------
    tuple[str, str, str]
        tuple[meneco_out, meneco_tsv, meneco_filtered]
        Path to the meneco json output, meneco tsv output and the filtered meneco tsv output corresponding to the
        meneco results of a step run (corresponding to the num)
    """
    meneco_out = os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{num}_meneco.json')
    meneco_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, f'{num}_meneco_out.tsv')
    meneco_filtered = os.path.join(OUTPUT, MENECO_D, FILTERED_D, f'{num}_meneco_out_filtered.tsv')
    return meneco_out, meneco_tsv, meneco_filtered


def run_meneco(network: str, output: str):
    """ Run Meneco tool on a SBML network and write the result in a JSON output.
    Use seeds, targets and bdd from Inputs directories.

    Parameters
    ----------
    network: str
        Network in SBML format to run the gapfilling toll Meneco on
    output: str
        Path to the JSON file to store meneco results
    """
    res = meneco.run_meneco(draftnet=network, seeds=SEEDS_ARTEFACTS[SBML_D], targets=TARGETS[SBML_D], repairnet=DB_SBML,
                            enumeration=False, json_output=True)

    json_output = json.dumps(res, indent=4)
    with open(output, "w") as outfile:
        outfile.write(json_output)


def meneco_json_to_tsv(output_json: str, output_tsv: str):
    """ Convert the JSON Meneco output to a TSV Meneco output thanks to padmet enhanced_meneco_output function.

    Parameters
    ----------
    output_json: str
        Path to Meneco output in JSON format as the input file.
    output_tsv: str
        Path to Meneco output in TSV format as the output file to create.
    """
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    os.system(f'padmet enhanced_meneco_output '
              f'--meneco={output_json} '
              f'--padmetRef={db} '
              f'--output={output_tsv} '
              f'--json '
              f'-v')


def add_rxn_to_nw(prev_nw: str, gap_filled_nw: str, rxn_to_add: str, rxn_to_exclude: Set[str] = None):
    """ Add reactions from Meneco in TSV format to a padmet network thanks to padmet manual_curation function.

    Parameters
    ----------
    prev_nw: str
        Network to add reaction in
    gap_filled_nw: str
        Output Network after adding reactions
    rxn_to_add: str
        TSV file with reaction to add to the network
    rxn_to_exclude: Set[str]
        Set of reactions to exclude (not added to network)
    """
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    if rxn_to_exclude is not None:
        with open(rxn_to_add, 'r') as f:
            lines = f.readlines()
        with open(rxn_to_add, 'w') as f:
            for line in lines:
                if line.split('\t')[0] not in rxn_to_exclude:
                    f.write(line)
    os.system(f'padmet manual_curation '
              f'--padmetSpec={prev_nw} '
              f'--data={rxn_to_add} '
              f'--padmetRef={db} '
              f'--output={gap_filled_nw} '
              f'--tool=MENECO '
              f'--category=GAP-FILLING '
              f'-v')


def extract_genes_from_blast() -> Dict[str, Set[str]]:
    """ Returns a dictionary associating for each reaction, all gene associated found with Blast hit


    Returns
    -------
    Dict[str, Set[str]]
        Dict[reaction, Set[genes_associated]]
        Associate for each reaction, all gene associated found with Blast hit
    """
    blast_res = os.path.join(OUTPUT, BLASTP_D, 'results', 'blast_results.tsv')
    dic_seq = dict()
    with open(blast_res, 'r') as blast_f:
        for line in blast_f:
            if not line.startswith('Reaction\t'):
                line = line.strip().split('\t')
                if line[0] not in dic_seq.keys():
                    dic_seq[line[0]] = set()
                dic_seq[line[0]].add(line[2])
    return dic_seq


def add_genes_tsv(tsv_file: str):
    """ Add the genes associated with each reaction in TSV file (meneco output filtered) for Blastp step.

    Parameters
    ----------
    tsv_file: str
        TSV file (meneco output filtered) to add the genes associated with each reaction in (for Blastp step)
    """
    logging.info(f'Adding genes linked to reactions in {tsv_file} file')

    dic_seq = extract_genes_from_blast()
    new_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, 'new_tsv.tsv')
    with open(tsv_file, 'r') as in_tsv, open(new_tsv, 'w') as out_tsv:
        for line in in_tsv:
            if line.startswith('idRef\t'):
                out_tsv.write(line.strip('\n'))
            else:
                line = line.strip().split('\t')
                genes = ' or '.join([f'({x})' for x in dic_seq[line[0]]])
                line.append(genes)
                out_tsv.write('\n' + '\t'.join(line))

    os.remove(tsv_file)
    os.rename(new_tsv, tsv_file)


def exists_not_producible_targets(output_json: str) -> bool:
    """ Returns if there still have not producible targets from meneco json output

    Parameters
    ----------
    output_json: str
        Meneco output in JSON format

    Returns
    -------
    bool
        True if there still have not producible targets from meneco json output
    """
    with open(output_json, 'r') as f:
        meneco_res = json.load(f)
    not_producible = meneco_res['Unproducible targets']
    return len(not_producible) > 0
