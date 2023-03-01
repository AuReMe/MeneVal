from environment import *
from typing import Tuple
import meneco
import json
import logging


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


def meneco_out_txt_to_tsv(output_txt: str, output_tsv: str):
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    os.system(f'padmet enhanced_meneco_output --meneco={output_txt} --padmetRef={db} --output={output_tsv} --json -v')


def add_rxn_to_nw(prev_nw: str, gap_filled_nw: str, rxn_to_add: str):
    db = get_file_from_ext(os.path.join(INPUT, DATABASE_D), PADMET_EXT)
    os.system(f'padmet manual_curation --padmetSpec={prev_nw} --data={rxn_to_add} --padmetRef={db} '
              f'--output={gap_filled_nw} --tool=MENECO --category=GAP-FILLING -v')
    print(f'{gap_filled_nw} created.')


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
