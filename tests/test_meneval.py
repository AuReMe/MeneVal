import os.path
import unittest
import shutil
from meneval.environment import *


class Test(unittest.TestCase):

    def setUp(self):
        os.system('meneval --init')
        self.move_files()
        os.system('meneval --check')

    # def tearDown(self):
    #     shutil.rmtree('Input')
    #     shutil.rmtree('Output')
    #     os.remove('meneco_validation.log')

    @staticmethod
    def move_files():
        shutil.rmtree('Input')
        source_dest_dir = {'Final_run/Input/': 'Input/'}
        source_dest_file = {'Final_run/Output/Networks/PADMET/0_base.padmet': 'Output/Networks/PADMET/0_base.padmet',
                            'Final_run/Output/Networks/SBML/0_base.sbml': 'Output/Networks/SBML/0_base.sbml'
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

        blastp_nw = get_nw_path(BLASTP)
        self.assertTrue(os.path.exists(blastp_nw[PADMET_D]))
        self.assertTrue(os.path.exists(blastp_nw[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 83)

    def test_enrich(self):
        groups = ['Group2', 'Group1', 'Group3']
        nb_lines = [83, 131, 179]

        for num in range(1, 4):
            group = groups[num-1]
            os.system(f'meneval --enrich {group}')

            self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, 'networks_validation.log')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, f'{group}_res_validation_networks.tsv')))

            self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, f'{num}_meneco_out_filtered.tsv')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{num}_meneco.json')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, f'{num}_meneco_out.tsv')))

            g_nw = get_nw_path(ENRICH, group)
            self.assertTrue(os.path.exists(g_nw[PADMET_D]))
            self.assertTrue(os.path.exists(g_nw[SBML_D]))

            with open('meneco_validation.log', 'r') as logfile:
                self.assertEqual(len(logfile.readlines()), nb_lines[num-1])


    def test_enrich_all(self):
        groups = ['Group2', 'Group1', 'Group3']
        group = 'ALL'
        os.system(f'meneval --enrich {group}')

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, 'networks_validation.log')))
        for g in groups:
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, f'{g}_res_validation_networks.tsv')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, f'{1}_meneco_out_filtered.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{1}_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, f'{1}_meneco_out.tsv')))

        g_nw = get_nw_path(ENRICH, group)
        self.assertTrue(os.path.exists(g_nw[PADMET_D]))
        self.assertTrue(os.path.exists(g_nw[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 83)


    def test_enrich_exclude(self):
        group = 'ALL'
        os.system(f'meneval --enrich {group}')
        os.system(f'meneval --fill')
        os.system(f'meneval --exclude')


    def test_fill(self):
        os.system('meneval --fill')

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, '1_meneco_out_filtered.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '1_meneco.json')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, '1_meneco_out.tsv')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_list.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_nb.tsv')))

        fill_nw = get_nw_path(FILL)
        self.assertTrue(os.path.exists(fill_nw[PADMET_D]))
        self.assertTrue(os.path.exists(fill_nw[SBML_D]))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 80)

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

        blastp_nw = get_nw_path(BLASTP)
        self.assertTrue(os.path.exists(blastp_nw[PADMET_D]))
        self.assertTrue(os.path.exists(blastp_nw[SBML_D]))

        # ENRICH
        groups = ['Group1', 'Group2', 'Group3']

        for num in range(2, 5):
            group = groups[num - 2]

            self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, 'networks_validation.log')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, ENRICH_D, group, 'res_validation_networks.tsv')))

            self.assertTrue(
                os.path.exists(os.path.join(OUTPUT, MENECO_D, FILTERED_D, f'{num}_meneco_out_filtered.tsv')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, f'{num}_meneco.json')))
            self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TSV_D, f'{num}_meneco_out.tsv')))

            g_nw = get_nw_path(ENRICH, group)
            self.assertTrue(os.path.exists(g_nw[PADMET_D]))
            self.assertTrue(os.path.exists(g_nw[SBML_D]))

        # FILL
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, TOOL_OUTPUTS_D, '5_meneco.json')))

        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_list.tsv')))
        self.assertTrue(os.path.exists(os.path.join(OUTPUT, MENECO_D, 'stat_nb.tsv')))

        with open('meneco_validation.log', 'r') as logfile:
            self.assertEqual(len(logfile.readlines()), 284)
