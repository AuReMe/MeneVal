import os.path
import unittest
import shutil

from aucomana.utils.reactions import Reactions

from meneval.validation_networks import *

# Parameters
RXN_LIST: List[str] = ['2-AMINOADIPATE-AMINOTRANSFERASE-RXN',
                       'HOMOCITRATE-SYNTHASE-RXN',
                       'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN',
                       'L-LYSINE-AMINOTRANSFERASE-RXN',
                       'ORNITHINE-CYCLODEAMINASE-RXN',
                       'RXN-13722',
                       'RXN-16756',
                       'RXN-19380',
                       'RXN-21797',
                       'RXN-21991',
                       'RXN-22438',
                       'RXN-5061',
                       'RXN-7970']
OUTPUT: str = 'Outputs_networks'
G1_REACTIONS_FILE: str = 'Final_run/Input/Enrichment/Group1/reactions.tsv'
G2_REACTIONS_FILE: str = 'Final_run/Input/Enrichment/Group2/reactions.tsv'
G3_REACTIONS_FILE: str = 'Final_run/Input/Enrichment/Group3/reactions.tsv'


class Test(unittest.TestCase):

    def setUp(self):
        os.mkdir(OUTPUT)

    def tearDown(self):
        shutil.rmtree(OUTPUT)

    def test_create_rxn_instance(self):
        rxn_group1 = create_rxn_instance(G1_REACTIONS_FILE, None, None)

        self.assertIsInstance(rxn_group1, Reactions)
        self.assertEqual(set(rxn_group1.species_list), {'ec042', 'HS'})
        self.assertEqual(rxn_group1.nb_reactions, 2361)

    def test_init_res_file(self):
        res_file = init_res_file('Group1', OUTPUT)
        self.assertTrue(os.path.exists(res_file))
        with open(res_file, 'r') as res_f:
            file_content = ['RXN\tNb Group1 presence\tGroup1 presence %\tGroup1 list']
            self.assertEqual(res_f.readlines(), file_content)

    def test_init_logger(self):
        init_logger(OUTPUT)
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'networks_validation.log')))

    def test_write_res(self):
        res_file = init_res_file('Group1', OUTPUT)

        # Write 1 result
        rxn = 'ORNITHINE-CYCLODEAMINASE-RXN'
        rxn_pres = ((2, 0.5), {'UTI89', 'HS'})
        write_res(rxn, rxn_pres, res_file)
        with open(res_file, 'r') as res_f:
            lines = res_f.readlines()
            exp_lines1 = ['RXN\tNb Group1 presence\tGroup1 presence %\tGroup1 list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tHS;UTI89']
            exp_lines2 = ['RXN\tNb Group1 presence\tGroup1 presence %\tGroup1 list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tUTI89;HS']
            self.assertTrue((lines == exp_lines1) or (lines == exp_lines2))

        # Write another result
        rxn = 'RXN-19380'
        rxn_pres = ((1, 0.25), {'IAI1'})
        write_res(rxn, rxn_pres, res_file)
        with open(res_file, 'r') as res_f:
            lines = res_f.readlines()
            exp_lines1 = ['RXN\tNb Group1 presence\tGroup1 presence %\tGroup1 list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tHS;UTI89\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            exp_lines2 = ['RXN\tNb Group1 presence\tGroup1 presence %\tGroup1 list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tUTI89;HS\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            self.assertTrue((lines == exp_lines1) or (lines == exp_lines2))

    def test_validation_networks(self):
        kept_rxn = validation_networks('Group2', OUTPUT, RXN_LIST, G2_REACTIONS_FILE)
        self.assertEqual(kept_rxn, {'RXN-19380'})
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'res_validation_networks.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'res_validation_networks.tsv')))
        with open(os.path.join(OUTPUT, 'res_validation_networks.tsv'), 'r') as res_f:
            lines = res_f.readlines()
            exp_lines = ['RXN\tNb Group2 presence\tGroup2 presence %\tGroup2 list\n',
                         'RXN-19380\t1\t50.0\tIAI1']
            self.assertEqual(lines, exp_lines)
