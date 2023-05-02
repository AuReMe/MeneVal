import os
import shutil
import unittest


class Test(unittest.TestCase):

    def setUp(self):
        os.system('meneval --init')

    def tearDown(self):
        shutil.rmtree('Input')
        shutil.rmtree('Output')
        os.remove('meneco_validation.log')

    def test_init(self):
        paths = ['Input',
                 'Input/AuCoMe',
                 'Input/DataBase',
                 'Input/Holobiont',
                 'Input/Species_seq',
                 'Input/Networks',
                 'Input/Seeds',
                 'Input/Targets',
                 'Output',
                 'Output/AuCoMe',
                 'Output/BlastP',
                 'Output/Holobiont',
                 'Output/Networks',
                 'Output/Networks/PADMET',
                 'Output/Networks/SBML',
                 'Output/Meneco',
                 'Output/Meneco/Filtered_TSV',
                 'Output/Meneco/Json_outputs',
                 'Output/Meneco/TSV']
        for file in paths:
            self.assertTrue(os.path.exists(file))

    @staticmethod
    def move_files():
        source_dest = {'Input_base/AuCoMe/group_template.tsv': 'Input/AuCoMe/group_template.tsv',
                       'Input_base/AuCoMe/reactions.tsv': 'Input/AuCoMe/reactions.tsv',
                       'Input_base/DataBase/metacyc_26.0_prot70.padmet': 'Input/DataBase/metacyc_26.0_prot70.padmet',
                       'Input_base/DataBase/proteins_seq_ids_reduced_70.fasta':
                           'Input/DataBase/proteins_seq_ids_reduced_70.fasta',
                       'Input_base/Networks/CFT073.padmet': 'Input/Networks/CFT073.padmet',
                       'Input_base/Seeds/artefacts.tsv': 'Input/Seeds/artefacts.tsv',
                       'Input_base/Seeds/seeds.tsv': 'Input/Seeds/seeds.tsv',
                       'Input_base/Species_seq/CFT073/CFT073.faa': 'Input/Species_seq/CFT073.faa',
                       'Input_base/Species_seq/CFT073/CFT073.fna': 'Input/Species_seq/CFT073.fna',
                       'Input_base/Targets/targets.tsv': 'Input/Targets/targets.tsv'
                       }
        for source, destination in source_dest.items():
            shutil.copy(source, destination)

    def test_check(self):
        self.move_files()
        os.system('meneval --check')

    def test_files_generator(self):
        self.move_files()
        files_to_be_create = ['Input/Seeds/seeds_medium.tsv',
                              'Input/Seeds/seeds_artefacts.tsv',
                              'Input/Seeds/seeds_artefacts.sbml',
                              'Input/Targets/temp_targets.tsv',
                              'Input/Targets/biomass.tsv',
                              'Input/Targets/targets.sbml',
                              'Input/DataBase/database.sbml',
                              'Output/Networks/PADMET/1_medium.padmet',
                              'Output/Networks/PADMET/1_base.padmet',
                              'Output/Networks/SBML/1_base.sbml']

        os.system('meneval --check')
        os.system('meneval --files')
        for file in files_to_be_create:
            self.assertTrue(os.path.exists(file))
