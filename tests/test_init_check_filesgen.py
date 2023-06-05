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
                 'Input/DataBase',
                 'Input/Enrichment',
                 'Input/Species_seq',
                 'Input/Networks',
                 'Input/Seeds',
                 'Input/Targets',
                 'Output',
                 'Output/BlastP',
                 'Output/Enrichment',
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
        source_dest = {'Final_run/Input/Enrichment/Group1/reactions.tsv': 'Input/Enrichment/Group1/reactions.tsv',
                       'Final_run/Input/Enrichment/Group2/reactions.tsv': 'Input/Enrichment/Group2/reactions.tsv',
                       'Final_run/Input/Enrichment/Group3/reactions.tsv': 'Input/Enrichment/Group3/reactions.tsv',
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
        dir_groups = ['Input/Enrichment/Group1', 'Input/Enrichment/Group2', 'Input/Enrichment/Group3']
        for dir_g in dir_groups:
            os.mkdir(dir_g)
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
                              'Output/Networks/PADMET/0_base.padmet',
                              'Output/Networks/SBML/0_base.sbml']

        os.system('meneval --check')
        os.system('meneval --files')
        for file in files_to_be_create:
            self.assertTrue(os.path.exists(file))
