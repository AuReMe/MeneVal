import os
import shutil
import unittest


class Test(unittest.TestCase):

    def setUp(self):
        os.system('python ../../meneval/meneval.py --init')

    def tearDown(self):
        self.addCleanup(shutil.rmtree, 'Input')
        self.addCleanup(shutil.rmtree, 'Output')
        self.addCleanup(os.remove, 'meneco_validation.log')

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
        source_dest = {'../Datas/AuCoMe/group_template.tsv': 'Input/AuCoMe/group_template.tsv',
                       '../Datas/AuCoMe/reactions.tsv': 'Input/AuCoMe/reactions.tsv',
                       '../Datas/DataBase/metacyc_26.0_prot70.padmet': 'Input/DataBase/metacyc_26.0_prot70.padmet',
                       '../Datas/DataBase/proteins_seq_ids_reduced_70.fasta':
                           'Input/DataBase/proteins_seq_ids_reduced_70.fasta',
                       '../Datas/Networks/CFT073.padmet': 'Input/Networks/CFT073.padmet',
                       '../Datas/Seeds/artefacts.tsv': 'Input/Seeds/artefacts.tsv',
                       '../Datas/Seeds/seeds.tsv': 'Input/Seeds/seeds.tsv',
                       '../Datas/Species_seq/CFT073/CFT073.faa': 'Input/Species_seq/CFT073.faa',
                       '../Datas/Species_seq/CFT073/CFT073.fna': 'Input/Species_seq/CFT073.fna',
                       '../Datas/Targets/targets.tsv': 'Input/Targets/targets.tsv'
                       }
        for source, destination in source_dest.items():
            shutil.copy(source, destination)

    def test_check(self):
        self.move_files()
        os.system('python ../../meneval/meneval.py --check')
