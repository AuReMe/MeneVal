import os.path
import unittest
import shutil
from meneval.environment import *


class Test(unittest.TestCase):

    def setUp(self):
        os.system('meneval --init')
        self.move_files()
        os.system('meneval --check')

    def tearDown(self):
        shutil.rmtree('Input')
        shutil.rmtree('Output')
        os.remove('meneco_validation.log')

    @staticmethod
    def move_files():
        shutil.rmtree('Input')
        source_dest_dir = {'Final_run/Input/': 'Input/'}
        source_dest_file = {'Final_run/Output/Networks/PADMET/1_base.padmet': 'Output/Networks/PADMET/1_base.padmet',
                            'Final_run/Output/Networks/PADMET/1_medium.padmet': 'Output/Networks/PADMET/1_medium.padmet',
                            'Final_run/Output/Networks/SBML/1_base.sbml': 'Output/Networks/SBML/1_base.sbml'
                            }
        for source, destination in source_dest_dir.items():
            shutil.copytree(source, destination)
        for source, destination in source_dest_file.items():
            shutil.copy(source, destination)

    def test_blastp(self):
        os.system('meneval --blastp')
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results', 'blast_results.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results', 'rxn_prot.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'sequences')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'blastp_validation.log')))
        seq_files = {'Q64602.fasta', 'Q2RJ84.fasta', 'Q07179.fasta', 'Q58667.fasta', 'Q9ZNE0.fasta', 'Q5JEW1.fasta',
                     'O59390.fasta', 'Q2RJ79.fasta', 'Q57926.fasta', 'Q8TW28.fasta', 'Q2RJ83.fasta', 'P70728.fasta',
                     'O59394.fasta', 'Q00852.fasta', 'P63510.fasta', 'Q44290.fasta', 'P54610.fasta', 'O94225.fasta',
                     'Q5SIJ1.fasta', 'O27668.fasta', 'Q01767.fasta', 'Q00853.fasta', 'O74298.fasta', 'Q8TLF1.fasta',
                     'J7SH14.fasta', 'Q8YMD9.fasta', 'O26917.fasta', 'Q57564.fasta', 'Q8TKQ6.fasta', 'Q58991.fasta',
                     'Q01181.fasta', 'Q58409.fasta', 'Q52070.fasta', 'P40976.fasta', 'P53090.fasta', 'Q59175.fasta',
                     'P07702.fasta', 'Q8TPT4.fasta', 'Q2RJ81.fasta', 'Q9ZND9.fasta', 'Q88H32.fasta', 'Q4J989.fasta',
                     'Q2RJ82.fasta', 'Q72LL6.fasta', 'Q2RJ80.fasta', 'P05342.fasta', 'P58350.fasta', 'P40495.fasta',
                     'NP_746228.1.fasta'}
        self.assertEqual(set(os.listdir(os.path.join(OUTPUT, BLASTP_D, 'sequences'))), seq_files)

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '1_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '1_meneco_out_filtered.tsv')))

        self.assertTrue(os.path.exists(BLASTP_GF_NW[PADMET_D]))
        self.assertTrue(os.path.exists(BLASTP_GF_NW[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 84)

    def test_aucome(self):
        os.system('meneval --aucome --group=group1')

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, AUCOME_D, 'networks_validation.log')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, AUCOME_D, 'res_validation_networks.tsv')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '3_meneco_out_filtered.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '3_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '3_meneco_out.tsv')))

        self.assertTrue(os.path.exists(AUCOME_GF_NW[PADMET_D]))
        self.assertTrue(os.path.exists(AUCOME_GF_NW[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 84)

    def test_fill(self):
        os.system('meneval --fill')

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '4_meneco_out_filtered.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '4_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '4_meneco_out.tsv')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_list.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_nb.tsv')))

        self.assertTrue(os.path.exists(FINAL_GF_NW[PADMET_D]))
        self.assertTrue(os.path.exists(FINAL_GF_NW[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 81)

    def test_workflow(self):
        os.system('meneval --workflow')

        # BLASTP
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results', 'blast_results.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'results', 'rxn_prot.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'sequences')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, BLASTP_D, 'blastp_validation.log')))
        seq_files = {'Q64602.fasta', 'Q2RJ84.fasta', 'Q07179.fasta', 'Q58667.fasta', 'Q9ZNE0.fasta', 'Q5JEW1.fasta',
                     'O59390.fasta', 'Q2RJ79.fasta', 'Q57926.fasta', 'Q8TW28.fasta', 'Q2RJ83.fasta', 'P70728.fasta',
                     'O59394.fasta', 'Q00852.fasta', 'P63510.fasta', 'Q44290.fasta', 'P54610.fasta', 'O94225.fasta',
                     'Q5SIJ1.fasta', 'O27668.fasta', 'Q01767.fasta', 'Q00853.fasta', 'O74298.fasta', 'Q8TLF1.fasta',
                     'J7SH14.fasta', 'Q8YMD9.fasta', 'O26917.fasta', 'Q57564.fasta', 'Q8TKQ6.fasta', 'Q58991.fasta',
                     'Q01181.fasta', 'Q58409.fasta', 'Q52070.fasta', 'P40976.fasta', 'P53090.fasta', 'Q59175.fasta',
                     'P07702.fasta', 'Q8TPT4.fasta', 'Q2RJ81.fasta', 'Q9ZND9.fasta', 'Q88H32.fasta', 'Q4J989.fasta',
                     'Q2RJ82.fasta', 'Q72LL6.fasta', 'Q2RJ80.fasta', 'P05342.fasta', 'P58350.fasta', 'P40495.fasta',
                     'NP_746228.1.fasta'}
        self.assertEqual(set(os.listdir(os.path.join(OUTPUT, BLASTP_D, 'sequences'))), seq_files)

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '1_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '1_meneco_out_filtered.tsv')))

        self.assertTrue(os.path.exists(BLASTP_GF_NW[PADMET_D]))
        self.assertTrue(os.path.exists(BLASTP_GF_NW[SBML_D]))

        # AUCOME
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, AUCOME_D, 'networks_validation.log')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, AUCOME_D, 'res_validation_networks.tsv')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '3_meneco_out_filtered.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '3_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '3_meneco_out.tsv')))

        self.assertTrue(os.path.exists(AUCOME_GF_NW[PADMET_D]))
        self.assertTrue(os.path.exists(AUCOME_GF_NW[SBML_D]))

        # FILL
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '4_meneco.json')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_list.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_nb.tsv')))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 188)
