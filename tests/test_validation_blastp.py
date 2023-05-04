import os.path
import unittest
import shutil
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
DB_PADMET: str = 'Final_run/Input/DataBase/metacyc_26.0_prot70.padmet'
PROT_FASTA_FILE: str = 'Final_run/Input/DataBase/proteins_seq_ids_reduced_70.fasta'
SPECIES_PROTEOME: str = 'Final_run/Input/Species_seq/CFT073.faa'
SPECIES_GENOME: str = 'Final_run/Input/Species_seq/CFT073.fna'

PADMET_SPEC: padmet.classes.padmetSpec = PadmetSpec(DB_PADMET)
PROT_FASTA = SeqIO.to_dict(SeqIO.parse(PROT_FASTA_FILE, 'fasta'))

UNIPROT_IDS = {'2-AMINOADIPATE-AMINOTRANSFERASE-RXN': {'UNIPROT:P58350', 'UNIPROT:Q64602', 'UNIPROT:P53090',
                                                       'PID:NP_746228.1', 'UNIPROT:Q72LL6'},
               'HOMOCITRATE-SYNTHASE-RXN': {'UNIPROT:Q00853', 'UNIPROT:P54610', 'UNIPROT:Q01181',
                                            'UNIPROT:Q8TKQ6', 'UNIPROT:Q4J989', 'UNIPROT:P05342',
                                            'UNIPROT:O59390', 'UNIPROT:P70728', 'UNIPROT:Q00852',
                                            'UNIPROT:O94225', 'UNIPROT:Q57926', 'UNIPROT:Q07179',
                                            'UNIPROT:Q44290', 'UNIPROT:Q8TW28', 'UNIPROT:Q52070'},
               'LYSINE--PYRUVATE-6-AMINOTRANSFERASE-RXN': set(),
               'L-LYSINE-AMINOTRANSFERASE-RXN': {'UNIPROT:Q01767', 'UNIPROT:P63510', 'UNIPROT:Q5JEW1'},
               'ORNITHINE-CYCLODEAMINASE-RXN': {'UNIPROT:Q59175', 'UNIPROT:Q88H32', 'UNIPROT:Q8YMD9',
                                                'UNIPROT:J7SH14'},
               'RXN-13722': {'UNIPROT:Q58667', 'UNIPROT:Q58409', 'UNIPROT:O27668', 'UNIPROT:Q9ZNE0',
                             'UNIPROT:Q9ZND9', 'UNIPROT:O26917', 'UNIPROT:Q8TLF1'},
               'RXN-16756': {'UNIPROT:P07702', 'UNIPROT:P40976', 'UNIPROT:O74298'},
               'RXN-19380': {'UNIPROT:Q57564', 'UNIPROT:Q8TPT4'},
               'RXN-21797': {'UNIPROT:Q2RJ79', 'UNIPROT:Q2RJ82', 'UNIPROT:Q2RJ84', 'UNIPROT:Q2RJ81',
                             'UNIPROT:Q2RJ83', 'UNIPROT:Q2RJ80'},
               'RXN-21991': set(),
               'RXN-22438': set(),
               'RXN-5061': set(),
               'RXN-7970': {'UNIPROT:P40495', 'UNIPROT:Q58991', 'UNIPROT:Q5SIJ1', 'UNIPROT:O59394'}}


class Test(unittest.TestCase):

    def setUp(self):
        os.mkdir(OUTPUT)

    def tearDown(self):
        shutil.rmtree(OUTPUT)

    def test_get_directories(self):
        expected_directories = ('Output_blastP/sequences',
                                'Output_blastP/results',
                                'Output_blastP/results/blast_results.tsv',
                                'Output_blastP/results/rxn_prot.tsv',
                                'Output_blastP/blastp_validation.log')
        directories = get_directories(OUTPUT)
        self.assertEqual(expected_directories, directories)

    def test_get_uniprot_ids_from_rxn(self):
        for rxn in RXN_LIST:
            uniprot_ids = get_uniprot_ids_from_rxn(rxn, PADMET_SPEC)
            self.assertEqual(uniprot_ids, UNIPROT_IDS[rxn])

    def test_create_dirs_and_init_result_file(self):
        seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(OUTPUT)
        create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file)
        dirs_to_create = [os.path.join(OUTPUT, 'results'),
                          os.path.join(OUTPUT, 'sequences'),
                          os.path.join(OUTPUT, 'results', 'blast_results.tsv')]
        for file in dirs_to_create:
            self.assertTrue(os.path.exists(file))
        with open(dirs_to_create[2], 'r') as blast_res:
            expected_line = ['Reaction\tUniprot ID\tSequence\tE value\tBit score\tIdentity (%)\tLength\tBlast method\n']
            lines = blast_res.readlines()
            self.assertEqual(expected_line, lines)

    def test_create_rxn_prot_tsv(self):
        seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(OUTPUT)
        create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file)
        create_rxn_prot_tsv(UNIPROT_IDS, rxn_prot_file)
        self.assertTrue(os.path.exists(rxn_prot_file))
        with open(rxn_prot_file, 'r') as r_prot:
            nb_lines = len(r_prot.readlines())
            self.assertEqual(nb_lines, len(RXN_LIST) + 1)

    def test_get_uniprot_seq(self):
        seq = {'UNIPROT:P58350': 'MTINATVKEAGFQPASRISSIGVSEILKIGARAAAMKREGKPVIILGAGEPDFDTPEHVKQAASDAIHRGETKYTALDGTPELKK'
                                 'AIREKFQRENGLAYELDEITVATGAKQILFNAMMASLDPGDEVIIPTPYWTSYSDIVHICEGKPVLIACDASSGFRLTAEKLEAA'
                                 'ITPRTRWVLLNSPSNPSGAAYSAADYRPLLEVLLRHPHVWLLVDDMYEHIVYDGFRFVTPAQLEPGLKNRTLTVNGVSKAYAMTG'
                                 'WRIGYAGGPRELIKAMAVVQSQATSCPSSISQAASVAALNGPQDFLKERTESFQRRRDLVVNGLNAIDGLDCRVPEGAFYTFSGC'
                                 'AGVLGKVTPSGKRIKTDTDFCAYLLEDAHVAVVPGSAFGLSPFFRISYATSEAELKEALERIAAACDRLS',
               'UNIPROT:P53090': 'MTLPESKDFSYLFSDETNARKPSPLKTCIHLFQDPNIIFLGGGLPLKDYFPWDNLSVDSPKPPFPQGIGAPIDEQNCIKYTVNKD'
                                 'YADKSANPSNDIPLSRALQYGFSAGQPELLNFIRDHTKIIHDLKYKDWDVLATAGNTNAWESTLRVFCNRGDVILVEAHSFSSSL'
                                 'ASAEAQGVITFPVPIDADGIIPEKLAKVMENWTPGAPKPKLLYTIPTGQNPTGTSIADHRKEAIYKIAQKYDFLIVEDEPYYFLQ'
                                 'MNPYIKDLKEREKAQSSPKQDHDEFLKSLANTFLSLDTEGRVIRMDSFSKVLAPGTRLGWITGSSKILKPYLSLHEMTIQAPAGF'
                                 'TQVLVNATLSRWGQKGYLDWLLGLRHEYTLKRDCAIDALYKYLPQSDAFVINPPIAGMFFTVNIDASVHPEFKTKYNSDPYQLEQ'
                                 'SLYHKVVERGVLVVPGSWFKSEGETEPPQPAESKEVSNPNIIFFRGTYAAVSPEKLTEGLKRLGDTLYEEFGISK',
               'UNIPROT:Q72LL6': 'MKPLSWSEAFGKGAGRIQASTIRELLKLTQRPGILSFAGGLPAPELFPKEEAAEAAARILREKGEVALQYSPTEGYAPLRAFVAE'
                                 'WIGVRPEEVLITTGSQQALDLVGKVFLDEGSPVLLEAPSYMGAIQAFRLQGPRFLTVPAGEEGPDLDALEEVLKRERPRFLYLIP'
                                 'SFQNPTGGLTPLPARKRLLQMVMERGLVVVEDDAYRELYFGEARLPSLFELAREAGYPGVIYLGSFSKVLSPGLRVAFAVAHPEA'
                                 'LQKLVQAKQGADLHTPMLNQMLVHELLKEGFSERLERVRRVYREKAQAMLHALDREVPKEVRYTRPKGGMFVWMELPKGLSAEGL'
                                 'FRRALEENVAFVPGGPFFANGGGENTLRLSYATLDREGIAEGVRRLGRALKGLLALV',
               'UNIPROT:Q64602': 'MNYSRFLTATSLARKTSPIRATVEIMSRAPKDIISLAPGSPNPKVFPFKSAVFTVENGSTIRFEGEMFQRALQYSSSYGIPELLS'
                                 'WLKQLQIKLHNPPTVNYSPNEGQMDLCITSGCQDGLCKVFEMLINPGDTVLVNEPLYSGALFAMKPLGCNFISVPSDDCGIIPEG'
                                 'LKKVLSQWKPEDSKDPTKRTPKFLYTIPNGNNPTGNSLTGDRKKEIYELARKYDFLIIEDDPYYFLQFTKPWEPTFLSMDVDGRV'
                                 'IRADSLSKVISSGLRVGFITGPKSLIQRIVLHTQISSLHPCTLSQLMISELLYQWGEEGFLAHVDRAIDFYKNQRDFILAAADKW'
                                 'LRGLAEWHVPKAGMFLWIKVNGISDAKKLIEEKAIEREILLVPGNSFFVDNSAPSSFFRASFSQVTPAQMDLVFQRLAQLIKDV'
                                 'S',
               'PID:NP_746228.1': 'MNQESISQSIAIVHPITLSHGRNAEVWDTDGKRYIDFVGGIGVLNLGHCNPAVVEAIQAQATRLTHYAFNAAPHGPYLALMEQL'
                                  'SQFVPVSYPLAGMLTNSGAEAAENALKVARGATGKRAIIAFDGGFHGRTLATLNLNGKVAPYKQRVGELPGPVYHLPYPSADTG'
                                  'VTCEQALKAMDRLFSVELAVEDVAAFIFEPVQGEGGFLALDPAFAQALRRFCDERGILIIIDEIQSGFGRTGQRFAFPRLGIEP'
                                  'DLLLLAKSIAGGMPLGAVVGRKELMAALPKGGLGGTYSGNPISCAAALASLAQMTDENLATWGERQEQAIVSRYERWKASGLSP'
                                  'YIGRLTGVGAMRGIEFANADGSPAPAQLAKVMEAARARGLLLMPSGKARHIIRLLAPLTIEAEVLEEGLDILEQCLAELN'
               }
        for uni_id in UNIPROT_IDS['2-AMINOADIPATE-AMINOTRANSFERASE-RXN']:
            self.assertEqual(get_uniprot_seq(uni_id, PROT_FASTA), seq[uni_id])

    def test_run_blastp(self):
        seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(OUTPUT)
        create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file)
        rxn = '2-AMINOADIPATE-AMINOTRANSFERASE-RXN'
        uni_id = 'UNIPROT:Q72LL6'
        seq = 'MKPLSWSEAFGKGAGRIQASTIRELLKLTQRPGILSFAGGLPAPELFPKEEAAEAAARILREKGEVALQYSPTEGYAPLRAFVAEWIGVRPEEVLITTGSQ' \
              'QALDLVGKVFLDEGSPVLLEAPSYMGAIQAFRLQGPRFLTVPAGEEGPDLDALEEVLKRERPRFLYLIPSFQNPTGGLTPLPARKRLLQMVMERGLVVVED' \
              'DAYRELYFGEARLPSLFELAREAGYPGVIYLGSFSKVLSPGLRVAFAVAHPEALQKLVQAKQGADLHTPMLNQMLVHELLKEGFSERLERVRRVYREKAQA' \
              'MLHALDREVPKEVRYTRPKGGMFVWMELPKGLSAEGLFRRALEENVAFVPGGPFFANGGGENTLRLSYATLDREGIAEGVRRLGRALKGLLALV'

        res = run_blastp(seq, uni_id, rxn, SPECIES_PROTEOME, seq_dir, blast_res_file)
        self.assertTrue(res)
        self.assertTrue(os.path.exists(os.path.join(seq_dir, 'Q72LL6.fasta')))

        exp_line = ['Reaction\tUniprot ID\tSequence\tE value\tBit score\tIdentity (%)\tLength\tBlast method\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tc1863\t1.56e-46\t163\t30.890\t382\tBlastp\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tc3224\t2.03e-24\t101\t27.200\t375\tBlastp\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tc2916\t3.79e-15\t73.6\t23.824\t340\tBlastp\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tc2548\t9.03e-12\t63.2\t27.778\t234\tBlastp\n']

        with open(blast_res_file, 'r') as res_file:
            self.assertEqual(res_file.readlines(), exp_line)

    def test_run_tblastn(self):
        seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(OUTPUT)
        create_dirs_and_init_result_file(seq_dir, res_dir, blast_res_file)
        rxn = '2-AMINOADIPATE-AMINOTRANSFERASE-RXN'
        uni_id = 'UNIPROT:Q72LL6'
        seq = 'MKPLSWSEAFGKGAGRIQASTIRELLKLTQRPGILSFAGGLPAPELFPKEEAAEAAARILREKGEVALQYSPTEGYAPLRAFVAEWIGVRPEEVLITTGSQ' \
              'QALDLVGKVFLDEGSPVLLEAPSYMGAIQAFRLQGPRFLTVPAGEEGPDLDALEEVLKRERPRFLYLIPSFQNPTGGLTPLPARKRLLQMVMERGLVVVED' \
              'DAYRELYFGEARLPSLFELAREAGYPGVIYLGSFSKVLSPGLRVAFAVAHPEALQKLVQAKQGADLHTPMLNQMLVHELLKEGFSERLERVRRVYREKAQA' \
              'MLHALDREVPKEVRYTRPKGGMFVWMELPKGLSAEGLFRRALEENVAFVPGGPFFANGGGENTLRLSYATLDREGIAEGVRRLGRALKGLLALV'

        res = run_tblastn(seq, uni_id, rxn, SPECIES_GENOME, seq_dir, blast_res_file)
        self.assertTrue(res)
        self.assertTrue(os.path.exists(os.path.join(seq_dir, 'Q72LL6.fasta')))

        exp_line = ['Reaction\tUniprot ID\tSequence\tE value\tBit score\tIdentity (%)\tLength\tBlast method\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tAE014075.1\t3.84e-38\t144\t32.222\t360\tTBlastN\n',
                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN\tQ72LL6\tAE014075.1\t1.64e-18\t85.5\t26.997\t363\tTBlastN\n']

        with open(blast_res_file, 'r') as res_file:
            self.assertEqual(res_file.readlines(), exp_line)

    def test_validation_blastp(self):
        kept_rxn = validation_blastp(RXN_LIST, OUTPUT, DB_PADMET, PROT_FASTA_FILE, SPECIES_PROTEOME, SPECIES_GENOME)
        print(kept_rxn)
        seq_dir, res_dir, blast_res_file, rxn_prot_file, log_file = get_directories(OUTPUT)
        seq_fasta = {'Q64602.fasta', 'Q2RJ84.fasta', 'Q07179.fasta', 'Q58667.fasta', 'Q9ZNE0.fasta', 'Q5JEW1.fasta',
                     'O59390.fasta', 'Q2RJ79.fasta', 'Q57926.fasta', 'Q8TW28.fasta', 'Q2RJ83.fasta', 'P70728.fasta',
                     'O59394.fasta', 'Q00852.fasta', 'P63510.fasta', 'Q44290.fasta', 'P54610.fasta', 'O94225.fasta',
                     'Q5SIJ1.fasta', 'O27668.fasta', 'Q01767.fasta', 'Q00853.fasta', 'O74298.fasta', 'Q8TLF1.fasta',
                     'J7SH14.fasta', 'Q8YMD9.fasta', 'O26917.fasta', 'Q57564.fasta', 'Q8TKQ6.fasta', 'Q58991.fasta',
                     'Q01181.fasta', 'Q58409.fasta', 'Q52070.fasta', 'P40976.fasta', 'P53090.fasta', 'Q59175.fasta',
                     'P07702.fasta', 'Q8TPT4.fasta', 'Q2RJ81.fasta', 'Q9ZND9.fasta', 'Q88H32.fasta', 'Q4J989.fasta',
                     'Q2RJ82.fasta', 'Q72LL6.fasta', 'Q2RJ80.fasta', 'P05342.fasta', 'P58350.fasta', 'P40495.fasta',
                     'NP_746228.1.fasta'}

        self.assertEqual(kept_rxn, {'L-LYSINE-AMINOTRANSFERASE-RXN', 'RXN-16756', 'RXN-13722', 'RXN-21797',
                                    '2-AMINOADIPATE-AMINOTRANSFERASE-RXN', 'HOMOCITRATE-SYNTHASE-RXN', 'RXN-7970'})
        self.assertEqual(set(os.listdir(seq_dir)), seq_fasta)
        self.assertTrue(os.path.exists(blast_res_file))
        self.assertTrue(os.path.exists(rxn_prot_file))
        self.assertTrue(os.path.exists(log_file))

        with open(blast_res_file, 'r') as blast_f, open(rxn_prot_file, 'r') as rxn_f, open(log_file, 'r') as log_f:
            self.assertEqual(len(blast_f.readlines()), 110)
            self.assertEqual(len(rxn_f.readlines()), 14)
            self.assertEqual(len(log_f.readlines()), 263)
