from os import path
from json import load

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align

# Fix AA compatibility

class Binder:
    def __init__(self, seq_rec: SeqRecord) -> None:
        seq = seq_rec._seq
        for k, v in seq_rec.__dict__.items():
            if not k.startswith('_'):
                self.__setattr__(k, v)

        self.phred = seq_rec._per_letter_annotations.get('phred_quality', [])
        self.mm_lmt = 6
        self.mode = self._select_mode(seq)
        self.seq = self._load_seq(seq)
        
        # Check reverse complement nt sequence if forward failed
        if self.seq is None and self.mode == 'nt':
            revcomp = Seq(seq).reverse_complement()
            self.seq = self._load_seq(revcomp)
        self.rfs = self._extract_rfs()

    @classmethod
    def seq(cls, seq: str):
        return cls(SeqRecord(seq))

    def _select_mode(self, seq):
        seq_set = set(seq)

        def _isnucleotide(seq_set) -> bool:
            nucleotides = set('ATGCU')
            return seq_set.issubset(nucleotides)

        def _isaminoacid(seq_set) -> bool:
            aminoacids = set('ARNDCQEGHILKMNFPOSUTWYVBZXJ*')
            return seq_set.issubset(aminoacids)

        if _isnucleotide(seq_set):
            return 'nt'
        elif _isaminoacid(seq_set):
            return 'aa'
        else:
            raise ValueError()

    def _load_rfs(self) -> dict:
        # TODO fix !!!
        def _translate_to_aa(rf_dict):
            for rf, params in rf_dict.items():
                for key, value in params.items():
                    if isinstance(value, str):
                        rf_dict[rf][key] = str(Seq(value).translate())
                    elif isinstance(value, int):
                        rf_dict[rf][key] = int(value / 3)
            return rf_dict
        
        file_dir = path.dirname(path.abspath(__file__))
        file_name = 'randomized_fragments.json'
        file_path = path.join(file_dir, file_name)
        with open(file_path, 'r') as file:
            rf_dict = load(file)

        rf_dict = rf_dict.get(self.lib)

        if self.mode == 'aa':
            return _translate_to_aa(rf_dict)
        return rf_dict

    def _load_seq(self, seq: Seq) -> Seq:

        # Remove non-letter characters except * (STOP)
        # def seq_cleanup(seq):
        #     chars = ascii_letters + '*'
        #     return ''.join([x for x in seq if x in chars])

        # seq = seq_cleanup(seq)

        # Empty string passed as seq
        if not seq:
            return None

        # TODO Extend tag lists/objects and create algorithm
        pelB = Seq('GCGGCCCAGCCGGCCATGGCG')
        his_tag = Seq('CACCATCACCACCATCAT')
        myc_tag = Seq('GAACAAAAACTCATCTCA')

        if self.mode == 'aa':
            pelB = pelB.translate()
            his_tag = his_tag.translate()
            myc_tag = myc_tag.translate()

        pelB_i = self._get_position(pelB, seq)
        his_tag_i = self._get_position(his_tag, seq)
        myc_tag_i = self._get_position(myc_tag, seq)

        end_i = max(his_tag_i, myc_tag_i)

        # Add length if sequence found
        if end_i > -1:
            end_i += len(his_tag)  # TODO myc_tag

        # Extract binder sequence
        seq_extract = seq[pelB_i:end_i]

        # If any sequence extracted
        if seq_extract:
            self.phred = self.phred[pelB_i:end_i] if self.phred else []
            return seq_extract
        else:
            return None

    def _get_position(self, target, seq=None) -> int:
        if seq is None:
            seq = self.seq

        aligner = Align.PairwiseAligner()
        # aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        aligner.query_internal_open_gap_score = -100
        aligner.query_internal_extend_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        # aligner.query_left_open_gap_score = 0
        # aligner.query_left_extend_gap_score = 0
        # aligner.query_right_open_gap_score = 0
        # aligner.query_right_extend_gap_score = 0
        aligner.mismatch_score = -0.5
        # aligner.match_score = 1

        alignments = aligner.align(str(seq), target)
        target_alignment = str(alignments[0]).split()[1]

        mm_condition = target_alignment.count('.') <= self.mm_lmt
        len_condition = len(target_alignment.strip('-')) == len(target)

        if len_condition and mm_condition:
            try:
                target_i = target_alignment.replace('.', '|').index('|')
                return target_i
            except ValueError:
                return -1
        return -1

    # TODO translate cannot be properly used | always auto translate?
    def _extract_rfs(self) -> list[Seq]:
        self.lib = self._lib_recognition()
        self.rf_dict = self._load_rfs()

        if self.seq is not None:

            # TODO load exact frame

            rf_extracted = []
            last_pos = (0, 0)

            # Iterate over randomized fragments
            for rf_data in self.rf_dict.values():
                rf = RF(rf_data, self.mode)
                extract_seq, *pos = rf.extract(self.seq)

                if self.mode == 'nt':
                    extract_seq = extract_seq.translate()

                if pos[0] > last_pos[1]:
                    rf_extracted.append(extract_seq)
                    last_pos = pos
                else:
                    rf_extracted.append(Seq('-'))

            return rf_extracted
            
        else:
            return [Seq('-')] * len(self.rf_dict.values())

    # def translate_rfs(self):
    #     if self.mode == 'nt':
    #         self.rfs = [rf.translate() for rf in self.rfs]

    # TODO refactor
    def _lib_recognition(self) -> str:

        lib_patterns = {
            'nt':
                [
                    {
                        'anchor': 'GCCCAGCCGGCCATG',
                        'i_start': 31,
                        'i_stop': 34,
                        'patterns':
                            {
                            'GCT': ['SH_VH3-23_Vk1-39',
                                    'PureLibra'],
                            'GCC': ['AI_VH3-23_Vk3-20',
                                    'AI_VH3-23_Vk1-39'],
                            'AGG': ['AI_VH1-69_Vk3-20',
                                    'AI_VH1-69_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'GCCCAGCCGGCCATG',
                        'i_start': 68,
                        'i_stop': 72,
                        'patterns':
                            {
                            'TGCA': ['SH_VH3-23_Vk1-39',
                                     'PureLibra'],
                            'CGCA': ['AI_VH3-23_Vk3-20',
                                     'AI_VH3-23_Vk1-39'],
                            'CAAG': ['AI_VH1-69_Vk3-20',
                                     'AI_VH1-69_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
                        'i_start': 25,
                        'i_stop': 29,
                        'patterns':
                            {
                                'GTGT': ['AI_VH3-23_Vk3-20',
                                         'AI_VH1-69_Vk3-20'],
                                'CAGA': ['SH_VH3-23_Vk1-39',
                                         'AI_VH3-23_Vk1-39',
                                         'AI_VH1-69_Vk1-39'],
                                'ATTG': ['PureLibra']
                            },
                    },
                    {
                        'anchor': 'GGTGGCGGTGGATCGGGCGGTGGTGG',
                        'i_start': 69,
                        'i_stop': 75,
                        'patterns':
                            {
                                'CCGTGT': ['SH_VH3-23_Vk1-39'],
                                'AAGAGC': ['AI_VH3-23_Vk3-20',
                                           'AI_VH1-69_Vk3-20'],
                                'CAGAGT': ['AI_VH3-23_Vk1-39',
                                           'AI_VH1-69_Vk1-39'],
                                'GGAAAG': ['PureLibra']
                            },
                    },
                ],
            'aa': 
                [
                    {
                        'anchor': 'AAQPAMA',
                        'i_start': 8,
                        'i_stop': 13,
                        'patterns':
                            {
                                'AEVKK': ['AI_VH1-69_Vk3-20',
                                          'AI_VH1-69_Vk1-39'],
                                'GGLVQ': ['AI_VH3-23_Vk3-20',
                                          'AI_VH3-23_Vk1-39',
                                          'SH_VH3-23_Vk1-39',
                                          'PureLibra'],
                            },
                    },
                    {
                        'anchor': 'AAQPAMA',
                        'i_start': 15,
                        'i_stop': 23,
                        'patterns':
                            {
                                'SSVKVSCK': ['AI_VH1-69_Vk3-20',
                                             'AI_VH1-69_Vk1-39'],
                                'GSLRLSCA': ['AI_VH3-23_Vk3-20',
                                             'AI_VH3-23_Vk1-39',
                                             'SH_VH3-23_Vk1-39',
                                             'PureLibra'],
                            },
                    },
                    {
                        'anchor': 'GGGGSGGGGSGGGGS',
                        'i_start': 0,
                        'i_stop': 4,
                        'patterns':
                            {
                                'EIVL': ['AI_VH3-23_Vk3-20',
                                        'AI_VH1-69_Vk3-20'],
                                'DIQM': ['AI_VH1-69_Vk1-39'
                                         'AI_VH3-23_Vk1-39',
                                         'SH_VH3-23_Vk1-39'],
                                'TEIV': ['PureLibra'],
                            },
                    },
                    {
                        'anchor': 'GGGGSGGGGSGGGGS',
                        'i_start': 12,
                        'i_stop': 22,
                        'patterns':
                            {
                                'SLSPGERATL': ['AI_VH3-23_Vk3-20',
                                               'AI_VH1-69_Vk3-20',
                                               'PureLibra'],
                                'ASVGDRVTIT': ['AI_VH1-69_Vk1-39'
                                               'AI_VH3-23_Vk1-39',
                                               'SH_VH3-23_Vk1-39'],
                            },
                    },
                    {
                        'anchor': 'WYQQKPG',
                        'i_start': 0,
                        'i_stop': 4,
                        'patterns':
                            {
                                'QAPR': ['AI_VH3-23_Vk3-20',
                                        'AI_VH1-69_Vk3-20',
                                        'PureLibra'],
                                'KAPK': ['AI_VH1-69_Vk1-39'
                                         'AI_VH3-23_Vk1-39',
                                         'SH_VH3-23_Vk1-39'],
                            },
                    },
                ],
            }
        

        lib_score = {
                'SH_VH3-23_Vk1-39': 0,
                'AI_VH3-23_Vk3-20': 0,
                'AI_VH3-23_Vk1-39': 0,
                'AI_VH1-69_Vk3-20': 0,
                'AI_VH1-69_Vk1-39': 0,
                'PureLibra': 0
        }

        for elem in lib_patterns[self.mode]:
            try:
                pos = self._get_position(elem['anchor']) + len(elem['anchor'])
                frag = self.seq[pos+elem['i_start']:pos+elem['i_stop']]
                rec_lib = elem['patterns'].get(frag)
                if rec_lib:
                    for match in rec_lib:
                        lib_score[match] += 1
            except TypeError as e:
                print(e)

        lib_sorted = sorted(list(lib_score.keys()),
                            key=lambda x: lib_score[x],
                            reverse=True)
        top_lib, second_lib, *rest = lib_sorted


        if lib_score[top_lib] > lib_score[second_lib]:
            return top_lib
        else:
            return 'default'

    def check_unique(self, file):
        # unique = True/False
        # return unique
        pass

        # [binder for binder in binder_list if binder.unique]


class RF:
    def __init__(self, params: dict, mode) -> None:
    
        self.mode = mode
        self.start = params['start']
        self.st_dis = len(self.start) + params['start_dis']
        self.end = params['end']
        self.end_dis = params['end_dis']
        self.limit = params['limit']
        self.end_pos = -1
        if mode == 'nt':
            self.mm_lmt = 7
        elif mode == 'aa':
            self.mm_lmt = 1

    def _get_position(self, target, seq) -> int:
        aligner = Align.PairwiseAligner()
        aligner.query_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        aligner.mismatch_score = -0.5

        target_alignments = aligner.align(str(seq), target)
        target_alignment = str(target_alignments[0]).split()

        mm_condition = target_alignment.count('.') <= self.mm_lmt
        len_condition = len(target_alignment[1].strip('-')) == len(target)

        if len_condition and mm_condition:
            try:
                target_i = target_alignment[1].replace('.', '|').index('|')
                return target_i
            except ValueError:
                return -1
        return -1

    # Extract RF sequence
    def extract(self, seq) -> tuple:
        self.st_pos = self._get_position(self.start, seq)
        self.end_pos = self._get_position(self.end, seq)

        # All positions found -> adjust distance
        if self._check_positions():
            self.rf_start = self.st_pos + self.st_dis
            self.rf_end = self.end_pos + self.end_dis
            seq = Seq(seq[self.rf_start:self.rf_end])

            if self.limit < len(seq):
                seq = Seq('-')

            return (seq, self.rf_start, self.rf_end)

        else:
            return (Seq('-'), 0, 0)

    def _check_positions(self) -> bool:
        conditions = [
            self.st_pos >= 0,
            self.end_pos > self.st_pos,
            self.st_pos + self.st_dis < self.end_pos + self.end_dis,
        ]
        return all(conditions)


if __name__ == '__main__':
    # seq = '''CAAAAACGTCTCAGGCAAGAAGGAGCACCGGCATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAGGAAGTGGAGCAAGACCGGGTTCGACTATTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCACGTGGAGGCGGTTCAAGCGGATGTGGCAGCGGCGGTGGCGGGTCGACAGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCACGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGAACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAATTGTATACGTTGTAGTTGACTTTTGGCCTGTGGACCAAGCTGGAGATCTATACGAGCGGCCGCAGAGGCCGCAGAACAAAAAACTCGTCTCAGAAGACGATCTGGGAAAACCCCGGTTGCCGACGGATTCTACTGTTTGAAACTTGTTTAGGAAAACCCTCATACTGAAAATTCATTTACCAAACCCCTGGAAAGAAAGACTAAACTTTATAATCGTTTCGCCTAAACTATCAAGGCTTTTCTATGGAAAAGCTAACAAGCCTTGGGGTTTGTACTGAGGGACAAAAACCCCTATTTTATGAAACAGGGGACCAAATTTGGGTTTTCCAACCCCCAAAAAATAAAGGGGGGGGGTTTCTCAAACGGGTGAAAATTTCTGAGAGGGGGCCTTTTTGAGGAAGGGGGCGGTATAAAAAAAACCCCCCCGTAGAAAAGGGGGAAACCTCAATATCCGCGGGAATAATTTTATATAACCCCTTTTCCGGGGGTTTTTTTTCTCTTGGGTCTTTTGAAAAAA'''

    #Pure Libra
    # seq = '''CCAAAACGTCTCAGGCAAGAAGGAGCACCGGCATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAGTGCGCGATGAAGCAGTCGATCGACTATTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCAGCGGCGGTGGCGGGTCGACGGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAGAAGGACCGGCGCATGTAGACTTTTGGCCAGGGGACCAAGCTGGAGATCAAACGAGCGGCCGCAGGGGCCGCAGAACAAAAACTCATCTCAGAAGAGGATCTGGGAGACGCGGGTGGCGGCGGTTCTACTGTTGAAAGTTGTTTAGCAAAACCTCATACAGAAAATTCATTTACTAACGTCTGGAAAGACGACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGCTACAGGCGTTGTGGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACATGGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAAGGGGGGTGGCTCTGAGGGTGGCGGTTCTGAAGGTGGGCGGTTCTGAAGGTGGCGATACAAACCTCCTGAGATAGGGATACACCCATTTCCGGGCAAACCAAAACGTCTCAGGCAAGAAGGAGCACCGGCATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAGTGCGCGATGAAGCAGTCGATCGACTATTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCAGCGGCGGTGGCGGGTCGACGGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAGAAGGACCGGCGCATGTAGACTTTTGGCCAGGGGACCAAGCTGGAGATCAAACGAGCGGCCGCAGGGGCCGCAGAACAAAAACTCATCTCAGAAGAGGATCTGGGAGACGCGGGTGGCGGCGGTTCTACTGTTGAAAGTTGTTTAGCAAAACCTCATACAGAAAATTCATTTACTAACGTCTGGAAAGACGACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTGTCTGTGGAATGCTACAGGCGTTGTGGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACATGGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAAGGGGGGTGGCTCTGAGGGTGGCGGTTCTGAAGGTGGGCGGTTCTGAAGGTGGCGATACAAACCTCCTGAGATAGGGATACACCCATTTCCGGGCAAA'''
    # seq = '''TTTGCCCGGAAATGGGTGTATCCCTATCTCAGGAGGTTTGTATCGCCACCTTCAGAACCGCCCACCTTCAGAACCGCCACCCTCAGAGCCACCCCCCTTCATTTTCAGGGATAGCAAGCCCAATAGGAACCCATGTACCGTAACACTGAGTTTCGTCACCAGTACAAACCACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGAGGTTTTGCTAAACAACTTTCAACAGTAGAACCGCCGCCACCCGCGTCTCCCAGATCCTCTTCTGAGATGAGTTTTTGTTCTGCGGCCCCTGCGGCCGCTCGTTTGATCTCCAGCTTGGTCCCCTGGCCAAAAGTCTACATGCGCCGGTCCTTCTGCTGACAGTAATACACTGCAAAATCTTCAGGCTCCAGTCTGCTGATGGTGAGAGTGAAGTCTGTCCCAGACCCACTGCCACTGAACCTGTCTGGGATGCCAGTGGCCCTGCTGGATGCACCATAGATGAGGAGCCTGGGAGCCTGGCCAGGTTTCTGCTGGTACCAGGCTAAGTAGCTGCTGCTAACACTCTGACTGGCCCTGCAGGAGAGGGTGGCTCTTTCCCCTGGAGACAAAGACAGGGTGCCTGGAGACTGCGTCAACACAATTTCCGTCGACCCGCCACCGCCGCTGCCACCTCCGCCTGAACCGCCTCCACCTGAGGACACGGTGACCAGGGTTCCCTGGCCCCAATAGTCGATCGACTGCTTCATCGCGCACTTCGCACAGTAATATACGGCCGTGTCCTCGGCTCTCAGGCTGTTCATTTGCAGATACAGCGTGTTCTTGGAATTGTCTCTGGAGATGGTGAACCGGCCCTTCACGGAGTCTGCGTAGTATGTGCTACCACCACTACCACTAATAGCTGACACCCACTCCAGCCCCTTCCCTGGAGCCTGGCGGACCCAGCTCATGGCATAGCTGCTAAAGGTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCCCCCAGGCTGTACCAAGCCTCCCCCAGACTCCAACAGCTGCACCTCGGCCATGGCCGGCTGGGCCGCGAGTAATAACAATCCAGCGGCTGCCGTAGGCAATAGGTATTTCATGCCGGTGCTCCTTCTTGCCTGAGACGTTTTGGTTTGCCCGGAAATGGGTGTATCCCTATCTCAGGAGGTTTGTATCGCCACCTTCAGAACCGCCCACCTTCAGAACCGCCACCCTCAGAGCCACCCCCCTTCATTTTCAGGGATAGCAAGCCCAATAGGAACCCATGTACCGTAACACTGAGTTTCGTCACCAGTACAAACCACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGAGGTTTTGCTAAACAACTTTCAACAGTAGAACCGCCGCCACCCGCGTCTCCCAGATCCTCTTCTGAGATGAGTTTTTGTTCTGCGGCCCCTGCGGCCGCTCGTTTGATCTCCAGCTTGGTCCCCTGGCCAAAAGTCTACATGCGCCGGTCCTTCTGCTGACAGTAATACACTGCAAAATCTTCAGGCTCCAGTCTGCTGATGGTGAGAGTGAAGTCTGTCCCAGACCCACTGCCACTGAACCTGTCTGGGATGCCAGTGGCCCTGCTGGATGCACCATAGATGAGGAGCCTGGGAGCCTGGCCAGGTTTCTGCTGGTACCAGGCTAAGTAGCTGCTGCTAACACTCTGACTGGCCCTGCAGGAGAGGGTGGCTCTTTCCCCTGGAGACAAAGACAGGGTGCCTGGAGACTGCGTCAACACAATTTCCGTCGACCCGCCACCGCCGCTGCCACCTCCGCCTGAACCGCCTCCACCTGAGGACACGGTGACCAGGGTTCCCTGGCCCCAATAGTCGATCGACTGCTTCATCGCGCACTTCGCACAGTAATATACGGCCGTGTCCTCGGCTCTCAGGCTGTTCATTTGCAGATACAGCGTGTTCTTGGAATTGTCTCTGGAGATGGTGAACCGGCCCTTCACGGAGTCTGCGTAGTATGTGCTACCACCACTACCACTAATAGCTGACACCCACTCCAGCCCCTTCCCTGGAGCCTGGCGGACCCAGCTCATGGCATAGCTGCTAAAGGTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCCCCCAGGCTGTACCAAGCCTCCCCCAGACTCCAACAGCTGCACCTCGGCCATGGCCGGCTGGGCCGCGAGTAATAACAATCCAGCGGCTGCCGTAGGCAATAGGTATTTCATGCCGGTGCTCCTTCTTGCCTGAGACGTTTTGG'''
    # seq = '''AAQPAMAEVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKCAMKQSIDYWGQGTLVTVSSGGGGSGGGGSGGGGSTEIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQKDRRM*TFGQGTKLEIKRAAAGAAEQKLISEEDLGDAGGGGSTVESCLAKPHTENSFTNVWKDDKTLDRYANYEGCLWNATGVVVCTGDETQCYGTWVPIGLAIPENEGGWL*GWRF*RWAVLKVAIQTS*DRDTPISGQTKTSQARRSTGMKYLLPTAAAGLLLLAAQPAMAEVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKCAMKQSIDYWGQGTLVTVSSGGGGSGGGGSGGGGSTEIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQKDRRM*TFGQGTKLEIKRAAAGAAEQKLISEEDLGDAGGGGSTVESCLAKPHTENSFTNVWKDDKTLDRYANYEGCLWNATGVVVCTGDETQCYGTWVPIGLAIPENEGGWL*GWRF*RWAVLKVAIQTS*DRDTPISGQ*TFGQGTKLEIKRAAAGAAEQKLIS'''

    #Pure Libra
    # seq = '''ATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCGGCCCAGCCGGCCATGGCCGAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTGTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGAAGTCGGATATTGTAGTAGTACCAGCTTTGACCCCTGGGGCCAGGGAACCCTGGTCACCGTGTCCTCAGGTGGAGGCGGTTCAGGCGGAGGTGGCAGCGGCGGTGGCGGGTCGACGGAAATTGTGTTGACGCAGTCTCCAGGCACCCTGTCTTTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCAGCAGGGCCACTGGCATCCCAGACAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACTCTCACCATCAGCAGACTGGAGCCTGAAGATTTTGCAGTGTATTACTGTCAGCAGTATGGTAGCTCACGGGGGTACTTTGGCCAGGGGACCAAGCTGGAGATCAAACGAGCGGCCGCAGGGGCCGCAGAACAAAAACTCATCTCAGAAGAGGATCTGGGAGACGCG'''
    # seq = '''TTTGCCCGGAAATGGGTGTATCCCTATCTCAGGAGGTTTGTATCGCCACCTTCAGAACCGCCCACCTTCAGAACCGCCACCCTCAGAGCCACCCCCCTTCATTTTCAGGGATAGCAAGCCCAATAGGAACCCATGTACCGTAACACTGAGTTTCGTCACCAGTACAAACCACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGAGGTTTTGCTAAACAACTTTCAACAGTAGAACCGCCGCCACCCGCGTCTCCCAGATCCTCTTCTGAGATGAGTTTTTGTTCTGCGGCCCCTGCGGCCGCTCGTTTGATCTCCAGCTTGGTCCCCTGGCCAAAAGTCTACATGCGCCGGTCCTTCTGCTGACAGTAATACACTGCAAAATCTTCAGGCTCCAGTCTGCTGATGGTGAGAGTGAAGTCTGTCCCAGACCCACTGCCACTGAACCTGTCTGGGATGCCAGTGGCCCTGCTGGATGCACCATAGATGAGGAGCCTGGGAGCCTGGCCAGGTTTCTGCTGGTACCAGGCTAAGTAGCTGCTGCTAACACTCTGACTGGCCCTGCAGGAGAGGGTGGCTCTTTCCCCTGGAGACAAAGACAGGGTGCCTGGAGACTGCGTCAACACAATTTCCGTCGACCCGCCACCGCCGCTGCCACCTCCGCCTGAACCGCCTCCACCTGAGGACACGGTGACCAGGGTTCCCTGGCCCCAATAGTCGATCGACTGCTTCATCGCGCACTTCGCACAGTAATATACGGCCGTGTCCTCGGCTCTCAGGCTGTTCATTTGCAGATACAGCGTGTTCTTGGAATTGTCTCTGGAGATGGTGAACCGGCCCTTCACGGAGTCTGCGTAGTATGTGCTACCACCACTACCACTAATAGCTGACACCCACTCCAGCCCCTTCCCTGGAGCCTGGCGGACCCAGCTCATGGCATAGCTGCTAAAGGTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCCCCCAGGCTGTACCAAGCCTCCCCCAGACTCCAACAGCTGCACCTCGGCCATGGCCGGCTGGGCCGCGAGTAATAACAATCCAGCGGCTGCCGTAGGCAATAGGTATTTCATGCCGGTGCTCCTTCTTGCCTGAGACGTTTTGGTTTGCCCGGAAATGGGTGTATCCCTATCTCAGGAGGTTTGTATCGCCACCTTCAGAACCGCCCACCTTCAGAACCGCCACCCTCAGAGCCACCCCCCTTCATTTTCAGGGATAGCAAGCCCAATAGGAACCCATGTACCGTAACACTGAGTTTCGTCACCAGTACAAACCACAACGCCTGTAGCATTCCACAGACAGCCCTCATAGTTAGCGTAACGATCTAAAGTTTTGTCGTCTTTCCAGACGTTAGTAAATGAATTTTCTGTATGAGGTTTTGCTAAACAACTTTCAACAGTAGAACCGCCGCCACCCGCGTCTCCCAGATCCTCTTCTGAGATGAGTTTTTGTTCTGCGGCCCCTGCGGCCGCTCGTTTGATCTCCAGCTTGGTCCCCTGGCCAAAAGTCTACATGCGCCGGTCCTTCTGCTGACAGTAATACACTGCAAAATCTTCAGGCTCCAGTCTGCTGATGGTGAGAGTGAAGTCTGTCCCAGACCCACTGCCACTGAACCTGTCTGGGATGCCAGTGGCCCTGCTGGATGCACCATAGATGAGGAGCCTGGGAGCCTGGCCAGGTTTCTGCTGGTACCAGGCTAAGTAGCTGCTGCTAACACTCTGACTGGCCCTGCAGGAGAGGGTGGCTCTTTCCCCTGGAGACAAAGACAGGGTGCCTGGAGACTGCGTCAACACAATTTCCGTCGACCCGCCACCGCCGCTGCCACCTCCGCCTGAACCGCCTCCACCTGAGGACACGGTGACCAGGGTTCCCTGGCCCCAATAGTCGATCGACTGCTTCATCGCGCACTTCGCACAGTAATATACGGCCGTGTCCTCGGCTCTCAGGCTGTTCATTTGCAGATACAGCGTGTTCTTGGAATTGTCTCTGGAGATGGTGAACCGGCCCTTCACGGAGTCTGCGTAGTATGTGCTACCACCACTACCACTAATAGCTGACACCCACTCCAGCCCCTTCCCTGGAGCCTGGCGGACCCAGCTCATGGCATAGCTGCTAAAGGTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCCCCCAGGCTGTACCAAGCCTCCCCCAGACTCCAACAGCTGCACCTCGGCCATGGCCGGCTGGGCCGCGAGTAATAACAATCCAGCGGCTGCCGTAGGCAATAGGTATTTCATGCCGGTGCTCCTTCTTGCCTGAGACGTTTTGG'''
    # seq = '''AAQPAMAEVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKEVGYCSSTSFDPWGQGTLVTVSSGGGGSGGGGSGGGGSTEIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSRGYFGQGTKLEIKRAAAGAAEQKLIS'''

    # seq = '''MKYLLPTAAAGLLLLAAQPAMAQVQLVQSGAEVKKPGSSVKVSCKASGDTFSTFDISWVRQAPGQGLEWMGGIVPMNYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCAREGGNWGAFDIWGQGTLVTVSSGGGGSGGGGSGGGGSEIVLTQSPGTLSLSPGERATLSCRASQSISNNYLAWYQQKPGQAPRLLIYDASSRATGIPDRFSGSGSGTDFTLTISSLEPEDFAVYYCQLYGSPPPRYTFGQGTKVEIKGPGGQHHHHHHGAYPYDVPDYA*'''
    # seq = ''
    # binder = Binder(seq)
    
    from pathlib import Path
    from Bio import SeqIO
    
    file = Path.cwd() / 'Desktop' / 'fastaq_test.fastq'
    
    with open(file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            # attribs = [attrib for attrib in dir(record) if not attrib.startswith('_')]
            # for attrib in attribs:
            #     print(attrib, eval(f'record.{attrib}'))
            # phred = record.letter_annotations['phred_quality']
            # print(phred)
            seq_rec = record
    

    # sec_rec = SeqRecord(seq)
    
    binder = Binder(seq_rec)
    
    # binder.load_seq_aa(seq)
    # print(binder.seq)
    # print(binder.seq)
    # print(binder.lib)
    print(binder.mode)
    print(binder.lib)
    # print(binder.seq)
    # binder.translate_rfs()
    print(binder.rfs)
    print(binder.features)
    print(binder.phred)
    
    if Seq('') in binder.rfs:
        print('OK')
