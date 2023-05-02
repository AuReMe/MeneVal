import os.path
import unittest
import shutil
from meneval.meneco_utils import *

# Parameters
MENECO_TSV = 'Final_run/Output/Meneco/TSV/1_meneco_out.tsv'
OUTPUT = 'Output_meneco_utils'


class Test(unittest.TestCase):

    def setUp(self):
        os.mkdir(OUTPUT)

    # def tearDown(self):
    #     shutil.rmtree('Input')
    #     shutil.rmtree('Output')
    #     os.remove('meneco_validation.log')

    def test_extract_rxn_from_meneco(self):
        rxn = extract_rxn_from_meneco(MENECO_TSV)
        exp_rxn = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN', 'HOMOCITRATE-SYNTHASE-RXN',
                   'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN', 'L-LYSINE-AMINOTRANSFERASE-RXN',
                   'ORNITHINE-CYCLODEAMINASE-RXN', 'RXN-13722', 'RXN-16756', 'RXN-19380', 'RXN-21797', 'RXN-21991',
                   'RXN-22438', 'RXN-5061', 'RXN-7970'}
        self.assertEqual(set(rxn), exp_rxn)

    def test_create_new_meneco_tsv(self):
        new_meneco = os.path.join(OUTPUT, 'new_meneco.tsv')
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

    def test_(self):
        pass
