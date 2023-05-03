import os.path
import unittest
import shutil
from meneval.meneco_utils import *
from meneval.environment import *

# Parameters
OUTPUT_M_UTILS = 'Output_meneco_utils'


class Test(unittest.TestCase):

    def setUp(self):
        os.mkdir(OUTPUT_M_UTILS)
        self.move_files()
        os.system('meneval --init')
        os.system('meneval --check')

    def tearDown(self):
        shutil.rmtree('Input')
        shutil.rmtree('Output')
        os.remove('meneco_validation.log')
        shutil.rmtree(OUTPUT_M_UTILS)

    @staticmethod
    def move_files():
        source_dest = {'Final_run/Input/': 'Input/',
                       'Final_run/Output/': 'Output/',
                       }
        for source, destination in source_dest.items():
            shutil.copytree(source, destination)

    def test_extract_rxn_from_meneco(self):
        meneco_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')
        rxn = extract_rxn_from_meneco(meneco_tsv)
        exp_rxn = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN', 'HOMOCITRATE-SYNTHASE-RXN',
                   'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN', 'L-LYSINE-AMINOTRANSFERASE-RXN',
                   'ORNITHINE-CYCLODEAMINASE-RXN', 'RXN-13722', 'RXN-16756', 'RXN-19380', 'RXN-21797', 'RXN-21991',
                   'RXN-22438', 'RXN-5061', 'RXN-7970'}
        self.assertEqual(set(rxn), exp_rxn)

    def test_create_new_meneco_tsv(self):
        meneco_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')
        new_meneco = os.path.join(OUTPUT_M_UTILS, 'new_meneco.tsv')
        kept_rxn = {'RXN-13722', 'RXN-16756', 'RXN-22438'}
        message = 'Adding for test'
        create_new_meneco_tsv(meneco_tsv, kept_rxn, new_meneco, message)
        self.assertTrue(os.path.exists(new_meneco))

        rxn = set(extract_rxn_from_meneco(new_meneco))
        self.assertEqual(rxn, kept_rxn)

        with open(new_meneco, 'r') as nm_f:
            for line in nm_f:
                line = line.split('\t')
                if line[6] != 'Comment':
                    self.assertEqual(line[6], message)

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

    def test_meneco_json_to_tsv(self):
        input_json = os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '1_meneco.json')
        output_tsv = os.path.join(OUTPUT_M_UTILS, 'output.tsv')
        exp_output_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')
        meneco_json_to_tsv(input_json, output_tsv)
        self.assertTrue(os.path.exists(output_tsv))
        with open(output_tsv, 'r') as o, open(exp_output_tsv, 'r') as exp_o:
            lines = o.readlines()
            exp_lines = exp_o.readlines()
            self.assertEqual(lines, exp_lines)

    def test_add_rxn_to_nw(self):
        prev_nw = BASE_NW[PADMET_D]
        new_nw = os.path.join(OUTPUT_M_UTILS, 'new_nw.padmet')
        rxn_tsv = os.path.join(OUTPUT, MENECO_D, FILTERED_D, '1_meneco_out_filtered.tsv')
        add_rxn_to_nw(prev_nw, new_nw, rxn_tsv)
        self.assertTrue(os.path.exists(new_nw))
        # TODO : Find a way to check all rxn was added

    def test_extract_genes_from_blast(self):
        genes_dict = extract_genes_from_blast()
        exp_genes_dict = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN': {'c1863', 'c3224', 'c3828', 'c2831', 'c2548', 'c0853',
                                                                  'c0189', 'c3210', 'c2148', 'c0688', 'c4134', 'c2916'},
                          'HOMOCITRATE-SYNTHASE-RXN': {'c0091'},
                          'L-LYSINE-AMINOTRANSFERASE-RXN': {'c3828', 'c0853', 'c0189', 'c3210', 'c2148', 'c4134'},
                          'RXN-13722': {'c0089', 'c1745', 'c0848', 'c0147', 'c0087'},
                          'RXN-16756': {'c2455', 'c0046', 'c2424', 'c2097', 'c0673', 'c2458', 'c2470', 'c2459', 'c2461',
                                        'c0454'},
                          'RXN-21797': {'c3974', 'c1025', 'c3456', 'c2995', 'c3464', 'c2679'},
                          'RXN-7970': {'c1517', 'c0090'}}

        self.assertEqual(genes_dict, exp_genes_dict)

    def test_add_genes_tsv(self):
        exp_tsv_file = os.path.join(OUTPUT, MENECO_D, FILTERED_D, '1_meneco_out_filtered.tsv')

        meneco_tsv = os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')
        kept_rxn = set(extract_rxn_from_meneco(exp_tsv_file))
        message = 'Adding for test'
        new_meneco = os.path.join(OUTPUT_M_UTILS, 'new_meneco.tsv')
        create_new_meneco_tsv(meneco_tsv, kept_rxn, new_meneco, message)

        add_genes_tsv(new_meneco)
        with open(new_meneco, 'r') as o, open(exp_tsv_file, 'r') as exp_o:
            lines = o.readlines()
            exp_lines = exp_o.readlines()
            for i in range(1, len(lines)):
                genes = set(lines[i].strip().split('\t')[7].split(' or '))
                exp_genes = set(exp_lines[i].strip().split('\t')[7].split(' or '))
                self.assertEqual(genes, exp_genes)

    def test_(self):
        pass










