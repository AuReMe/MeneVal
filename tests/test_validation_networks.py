import os.path
import unittest
import shutil

import aucomana.utils.reactions

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
NAME_SPECIES: str = 'group1_name'
OUTPUT: str = 'Outputs_networks'
REACTIONS_FILE: str = 'Final_run/Input/AuCoMe/reactions.tsv'
GROUP_FILE: str = 'Final_run/Input/AuCoMe/group_template.tsv'
GROUP: str = 'group1'

GROUP1 = {'UTI89', 'IAI1', 'CFT073', 'HS'}
ALL_SP = {'HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82'}


class Test(unittest.TestCase):

    def setUp(self):
        os.mkdir(OUTPUT)

    def tearDown(self):
        shutil.rmtree(OUTPUT)

    def test_create_rxn_instance(self):
        rxn_group = create_rxn_instance(REACTIONS_FILE, GROUP_FILE, GROUP)
        rxn_no_group = create_rxn_instance(REACTIONS_FILE, None, None)
        rxn_group_no_group = create_rxn_instance(REACTIONS_FILE, GROUP_FILE, None)
        rxn_group_no_group_file = create_rxn_instance(REACTIONS_FILE, None, GROUP)

        self.assertIsInstance(rxn_group, aucomana.utils.reactions.Reactions)
        self.assertIsInstance(rxn_no_group, aucomana.utils.reactions.Reactions)
        self.assertIsInstance(rxn_group_no_group, aucomana.utils.reactions.Reactions)
        self.assertIsInstance(rxn_group_no_group_file, aucomana.utils.reactions.Reactions)

        self.assertEqual(set(rxn_group.species_list), GROUP1)
        self.assertEqual(set(rxn_no_group.species_list), ALL_SP)
        self.assertEqual(set(rxn_group_no_group.species_list), ALL_SP)
        self.assertEqual(set(rxn_group_no_group_file.species_list), ALL_SP)

    def test_init_res_file(self):
        res_file = init_res_file(NAME_SPECIES, OUTPUT)
        self.assertTrue(os.path.exists(res_file))
        with open(res_file, 'r') as res_f:
            file_content = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list']
            self.assertEqual(res_f.readlines(), file_content)

    def test_init_logger(self):
        init_logger(OUTPUT)
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'networks_validation.log')))

    def test_write_res(self):
        res_file = init_res_file(NAME_SPECIES, OUTPUT)

        # Write 1 result
        rxn = 'ORNITHINE-CYCLODEAMINASE-RXN'
        rxn_pres = ((2, 0.5), {'UTI89', 'HS'})
        write_res(rxn, rxn_pres, res_file)
        with open(res_file, 'r') as res_f:
            lines = res_f.readlines()
            exp_lines1 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tHS;UTI89']
            exp_lines2 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tUTI89;HS']
            self.assertTrue((lines == exp_lines1) or (lines == exp_lines2))

        # Write another result
        rxn = 'RXN-19380'
        rxn_pres = ((1, 0.25), {'IAI1'})
        write_res(rxn, rxn_pres, res_file)
        with open(res_file, 'r') as res_f:
            lines = res_f.readlines()
            exp_lines1 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tHS;UTI89\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            exp_lines2 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tUTI89;HS\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            self.assertTrue((lines == exp_lines1) or (lines == exp_lines2))

    def test_validation_networks(self):
        kept_rxn = validation_networks(NAME_SPECIES, OUTPUT, RXN_LIST, REACTIONS_FILE, GROUP_FILE, GROUP)
        self.assertEqual(kept_rxn, {'ORNITHINE-CYCLODEAMINASE-RXN', 'RXN-19380'})
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'res_validation_networks.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, 'res_validation_networks.tsv')))
        with open(os.path.join(OUTPUT, 'res_validation_networks.tsv'), 'r') as res_f:
            lines = res_f.readlines()
            exp_lines1 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tHS;UTI89\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            exp_lines2 = ['RXN\tNb group1_name presence\tgroup1_name presence %\tgroup1_name list\n',
                          'ORNITHINE-CYCLODEAMINASE-RXN\t2\t50.0\tUTI89;HS\n',
                          'RXN-19380\t1\t25.0\tIAI1']
            self.assertTrue((lines == exp_lines1) or (lines == exp_lines2))
