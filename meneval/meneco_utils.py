import meneco
import json
import logging
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


def create_new_meneco_tsv(meneco_tsv: str, kept_rxn: Set[str], output_meneco_tsv: str, message: str):
    """Create a new Meneco tsv output keeping only the reactions that led to a blast match.

    Parameters
    ----------
    meneco_tsv: str
        Meneco results in TSV format (output from enhanced_meneco_output in padmet)
    kept_rxn : Set[str]
        set of reactions to keep in the Meneco output.
    output_meneco_tsv: str
        Path to the new meneco tsv output
    message: str
        Comment for adding the reaction
    """
    logging.info('\nCreation of filtered Meneco tsv output file')
    with open(meneco_tsv, 'r') as in_meneco, open(output_meneco_tsv, 'w') as out_meneco:
        for line in in_meneco:
            if line.startswith('idRef'):
                out_meneco.write(line)
            elif line.split('\t')[0] in kept_rxn:
                line = line.split('\t')
                line[-2] = message
                out_meneco.write('\t'.join(line))
    logging.info(f'Meneco tsv output file saved in : {output_meneco_tsv}')


def get_meneco_files(num: int) -> Tuple[str, str, str]:
    meneco_out = os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{num}_meneco.json')
    meneco_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, f'{num}_meneco_out.tsv')
    meneco_filtered = os.path.join(OUTPUT, MENECO_D, FILTERED_D, f'{num}_meneco_out_filtered.tsv')
    return meneco_out, meneco_tsv, meneco_filtered


def run_meneco(network: str, output: str):
    res = meneco.run_meneco(draftnet=network, seeds=SEEDS_ARTEFACTS[SBML_D], targets=TARGETS[SBML_D], repairnet=DB_SBML,
                            enumeration=False, json_output=True)

    json_output = json.dumps(res, indent=4)
    with open(output, "w") as outfile:
        outfile.write(json_output)


def meneco_json_to_tsv(output_json: str, output_tsv: str):
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    os.system(f'padmet enhanced_meneco_output '
              f'--meneco={output_json} '
              f'--padmetRef={db} '
              f'--output={output_tsv} '
              f'--json '
              f'-v')


def add_rxn_to_nw(prev_nw: str, gap_filled_nw: str, rxn_to_add: str):
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    os.system(f'padmet manual_curation --padmetSpec={prev_nw} --data={rxn_to_add} --padmetRef={db} '
              f'--output={gap_filled_nw} --tool=MENECO --category=GAP-FILLING -v')


def extract_genes_from_blast():
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


def add_genes_tsv(tsv_file):
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


def exists_not_producible_targets(output_json):
    with open(output_json, 'r') as f:
        meneco_res = json.load(f)
    not_producible = meneco_res['Unproducible targets']
    return len(not_producible) > 0
