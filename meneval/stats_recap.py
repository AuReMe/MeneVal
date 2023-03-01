import json
from padmet.utils.connection.sbml_to_curation_form import convert_from_coded_id
from meneval.environment import *


HEADER_NB = ['Nb not producible', 'Nb not constructable', 'Nb constructible', 'Nb essential', 'Nb minimal',
             'Nb intersection', 'Nb union']

HEADER_L = ['Not producible', 'Not constructable', 'Constructible', 'Essential', 'Minimal', 'Intersection', 'Union']

ROWS = ['', 'Initial', 'After BlastP', 'After Holobiont', 'After Aucome']

FILE_NB = os.path.join(OUTPUT, MENECO_D, 'stat_nb.tsv')
FILE_L = os.path.join(OUTPUT, MENECO_D, 'stat_list.tsv')


def make_meneco_stats():
    with open(FILE_NB, 'w') as fnb, open(FILE_L, 'w') as fl:
        fnb.write('\t'.join([ROWS[0]] + HEADER_NB))
        fl.write('\t'.join([ROWS[0]] + HEADER_L))

        for i in range(1, 5):
            output = os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{i}_meneco.json')

            with open(output, 'r') as o:
                res = json.load(o)

                essential = res['Essential reactions']
                reactions_set = set()
                for l_rxn in essential.values():
                    for rxn in l_rxn:
                        reactions_set.add(rxn)
                essential_rxn = list(reactions_set)

                stat = [
                    [convert_from_coded_id(x)[0] for x in res['Unproducible targets']],
                    [convert_from_coded_id(x)[0] for x in res['Unreconstructable targets']],
                    [convert_from_coded_id(x)[0] for x in res['Reconstructable targets']],
                    [convert_from_coded_id(x)[0] for x in essential_rxn],
                    [convert_from_coded_id(x)[0] for x in res['One minimal completion']],
                    [convert_from_coded_id(x)[0] for x in res['Intersection of cardinality minimal completions']],
                    [convert_from_coded_id(x)[0] for x in res['Union of cardinality minimal completions']]
                ]

                fnb.write('\n' + ROWS[i] + '\t' + '\t'.join([str(len(x)) for x in stat]))
                fl.write('\n' + ROWS[i] + '\t' + '\t'.join([', '.join(x) for x in stat]))
