import unittest
from meneval.validation_BlastP import *

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
OUTPUT: str = 'Output_blastP'
DB_PADMET: str = 'Input_files_generated/DataBase/metacyc_26.0_prot70.padmet'
PROT_FASTA: str = 'Input_files_generated/DataBase/proteins_seq_ids_reduced_70.fasta'
SPECIES_PROTEOME: str = 'Input_files_generated/Species_seq/CFT073.faa'
SPECIES_GENOME: str = 'Input_files_generated/Species_seq/CFT073.fna'


class Test(unittest.TestCase):

    def test_get_directories(self):
        # os.mkdir(OUTPUT)
        os.system('meneval --blastp')
