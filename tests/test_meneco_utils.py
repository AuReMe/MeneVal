import os.path
import unittest
import shutil
from meneval.meneco_utils import *

# Parameters
MENECO_TSV = 'Final_run/Output/Meneco/TSV/1_meneco_out.tsv'
OUTPUT = 'Output_meneco_utils'


class Test(unittest.TestCase):

    def setUp(self):
        os.system('meneval --init')

    # def tearDown(self):
    #     shutil.rmtree('Input')
    #     shutil.rmtree('Output')
    #     os.remove('meneco_validation.log')

    @staticmethod
    def move_files():
        source_dest = {'Final_run/Input/AuCoMe/group_template.tsv': 'Input/AuCoMe/group_template.tsv',
                       'Final_run/Input/AuCoMe/reactions.tsv': 'Input/AuCoMe/reactions.tsv',
                       'Final_run/Input/DataBase/metacyc_26.0_prot70.padmet':
                           'Input/DataBase/metacyc_26.0_prot70.padmet',
                       'Final_run/Input/DataBase/proteins_seq_ids_reduced_70.fasta':
                           'Input/DataBase/proteins_seq_ids_reduced_70.fasta',
                       'Final_run/Input/Networks/CFT073.padmet': 'Input/Networks/CFT073.padmet',
                       'Final_run/Input/Seeds/artefacts.tsv': 'Input/Seeds/artefacts.tsv',
                       'Final_run/Input/Seeds/seeds.tsv': 'Input/Seeds/seeds.tsv',
                       'Final_run/Input/Species_seq/CFT073.faa': 'Input/Species_seq/CFT073.faa',
                       'Final_run/Input/Species_seq/CFT073.fna': 'Input/Species_seq/CFT073.fna',
                       'Final_run/Input/Targets/targets.tsv': 'Input/Targets/targets.tsv'
                       }
        for source, destination in source_dest.items():
            shutil.copy(source, destination)

    def test_extract_rxn_from_meneco(self):
        rxn = extract_rxn_from_meneco(MENECO_TSV)
        exp_rxn = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN', 'HOMOCITRATE-SYNTHASE-RXN',
                   'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN', 'L-LYSINE-AMINOTRANSFERASE-RXN',
                   'ORNITHINE-CYCLODEAMINASE-RXN', 'RXN-13722', 'RXN-16756', 'RXN-19380', 'RXN-21797', 'RXN-21991',
                   'RXN-22438', 'RXN-5061', 'RXN-7970'}
        self.assertEqual(set(rxn), exp_rxn)

    def test_create_new_meneco_tsv(self):
        os.mkdir(OUTPUT)
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
        shutil.rmtree(OUTPUT)

    def test_get_meneco_files(self):
        pass
