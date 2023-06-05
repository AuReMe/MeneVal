import csv
import json

from padmet.utils.sbmlPlugin import convert_from_coded_id

BASE = 'base.json'
FILLED = 'filled.json'
TARGETS = 'Input/Targets/targets.tsv'
OUTPUT = 'targets_scope.tsv'
MENECO = 'Output/Meneco/Json_outputs/4_meneco.json'


def check_prod():
    with open(BASE, 'r') as b, open(FILLED, 'r') as f:
        base_d = json.load(b)
        filled_d = json.load(f)

    for k, v in base_d.items():
        new_v = set()
        for met in v:
            new_met = convert_from_coded_id(met)[0]
            new_v.add(new_met)
        base_d[k] = new_v

    for k, v in filled_d.items():
        new_v = set()
        for met in v:
            new_met = convert_from_coded_id(met)[0]
            new_v.add(new_met)
        filled_d[k] = new_v

    with open(TARGETS, 'r') as it, open(OUTPUT, 'w') as ot:
        ot.write('\t'.join(['Metabolite', 'Scope Base', 'Scope Gapfilling']))
        it = csv.reader(it, delimiter='\t')
        for row in it:
            met = row[1]
            if met in base_d['producible_target']:
                base_prod = '1'
            else:
                base_prod = '0'
            if met in filled_d['producible_target']:
                filled_prod = '1'
            else:
                filled_prod = '0'
            ot.write('\n' + '\t'.join([met, base_prod, filled_prod]))


def extract_rxn():
    with open(MENECO, 'r') as m:
        res = json.load(m)['Essential reactions']
    for r, m_l in res.items():
        r = convert_from_coded_id(r)[0]
        print()
        print(r)
        for m in m_l:
            m = convert_from_coded_id(m)[0]
            print(m)



extract_rxn()