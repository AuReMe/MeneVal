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
OUTPUT: str = 'Ouputs_networks'
REACTIONS_FILE: str = 'Files_generated/Input/AuCoMe/reactions.tsv'
GROUP_FILE: str = 'Files_generated/Input/AuCoMe/group_template.tsv'
GROUP: str = 'group1'

GROUP1 = {'UTI89', 'IAI1', 'CFT073', 'HS'}
ALL_SP = {'HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82'}


class Test(unittest.TestCase):

    # def setUp(self):
    #     os.mkdir(OUTPUT)
    #
    # def tearDown(self):
    #     shutil.rmtree(OUTPUT)

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

    def test_(self):
        pass
