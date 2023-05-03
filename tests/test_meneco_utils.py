import os.path
import unittest
import shutil
from meneval.meneco_utils import *
from meneval.environment import *

# Parameters
MENECO_TSV = 'Final_run/Output/Meneco/TSV/1_meneco_out.tsv'
OUTPUT_M_UTILS = 'Output_meneco_utils'


class Test(unittest.TestCase):

    def setUp(self):
        self.move_files()
        os.system('meneval --init')
        os.system('meneval --check')

    def tearDown(self):
        shutil.rmtree('Input')
        shutil.rmtree('Output')
        os.remove('meneco_validation.log')

    @staticmethod
    def move_files():
        source_dest = {'Final_run/Input/': 'Input/',
                       'Final_run/Output/': 'Output/',
                       }
        for source, destination in source_dest.items():
            shutil.copytree(source, destination)

    def test_extract_rxn_from_meneco(self):
        rxn = extract_rxn_from_meneco(MENECO_TSV)
        exp_rxn = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN', 'HOMOCITRATE-SYNTHASE-RXN',
                   'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN', 'L-LYSINE-AMINOTRANSFERASE-RXN',
                   'ORNITHINE-CYCLODEAMINASE-RXN', 'RXN-13722', 'RXN-16756', 'RXN-19380', 'RXN-21797', 'RXN-21991',
                   'RXN-22438', 'RXN-5061', 'RXN-7970'}
        self.assertEqual(set(rxn), exp_rxn)

    def test_create_new_meneco_tsv(self):
        os.mkdir(OUTPUT_M_UTILS)
        new_meneco = os.path.join(OUTPUT_M_UTILS, 'new_meneco.tsv')
        kept_rxn = {'RXN-13722', 'RXN-16756', 'RXN-22438'}
        message = 'Adding for test'
        create_new_meneco_tsv(MENECO_TSV, kept_rxn, new_meneco, message)
        self.assertTrue(os.path.exists(new_meneco))

        rxn = set(extract_rxn_from_meneco(new_meneco))
        self.assertEqual(rxn, kept_rxn)

        with open(new_meneco, 'r') as nm_f:
            for line in nm_f:
                line = line.split('\t')
                if line[6] != 'Comment':
                    self.assertEqual(line[6], message)
        shutil.rmtree(OUTPUT_M_UTILS)

    def test_get_meneco_files(self):
        exp_num1 = ('Output/Meneco/Json_outputs/1_meneco.json',
                    'Output/Meneco/TSV/1_meneco_out.tsv',
                    'Output/Meneco/Filtered_TSV/1_meneco_out_filtered.tsv')
        exp_num4 = ('Output/Meneco/Json_outputs/4_meneco.json',
                    'Output/Meneco/TSV/4_meneco_out.tsv',
                    'Output/Meneco/Filtered_TSV/4_meneco_out_filtered.tsv')
        self.assertEqual(get_meneco_files(1), exp_num1)
        self.assertEqual(get_meneco_files(4), exp_num4)

    def test_run_meneco(self):
        os.mkdir(OUTPUT_M_UTILS)
        network = BASE_NW[SBML_D]
        output = os.path.join(OUTPUT_M_UTILS, 'output.json')
        run_meneco(network, output)
        self.assertTrue(os.path.exists(output))

        exp_output = os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '1_meneco.json')
        with open(output, 'r') as o, open(exp_output, 'r') as exp_o:
            res = json.load(o)
            exp_res = json.load(exp_o)
        self.assertEqual(res['Draft network file'], exp_res['Draft network file'])
        self.assertEqual(res['Seeds file'], exp_res['Seeds file'])
        self.assertEqual(res['Targets file'], exp_res['Targets file'])
        self.assertEqual(res['Repair db file'], exp_res['Repair db file'])
        self.assertEqual(set(res['Unproducible targets']), set(exp_res['Unproducible targets']))
        self.assertEqual(set(res['Unreconstructable targets']), set(exp_res['Unreconstructable targets']))
        self.assertEqual(set(res['Reconstructable targets']), set(exp_res['Reconstructable targets']))
        self.assertEqual(set(res['Intersection of cardinality minimal completions']),
                         set(exp_res['Intersection of cardinality minimal completions']))
        self.assertEqual(set(res['Union of cardinality minimal completions']),
                         set(exp_res['Union of cardinality minimal completions']))

        shutil.rmtree(OUTPUT_M_UTILS)

    def test_(self):
        pass








