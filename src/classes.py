from os import path
from json import load
from collections import Counter, defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import Align
import re
from typing import NamedTuple

from src.data import lib_frames, supE_codon_table, seq_liab_db
# from data import lib_frames, supE_codon_table


class Result(NamedTuple):
    seq_id: str
    rfs: list
    lib: str
    quality: tuple
    frame_muts: list
    stop_muts: list
    seq_aa: str


class Binder:
    def __init__(self, seq_rec: SeqRecord) -> None:

        self.mode = self._select_mode(seq_rec)
        if self.mode == 'nt':
            self.mm_lmt = 7
        elif self.mode == 'aa':
            self.mm_lmt = 1
        self.hq_lmt = 40  #

        if isinstance(seq_rec._seq, str):
            seq_rec._seq = Seq(seq_rec._seq)

        for k, v in seq_rec.__dict__.items():
            if not k.startswith('_'):
                self.__setattr__(k, v)

        seq_extract = self._extract_seq(seq_rec)

        # Check reverse complement nt sequence if forward failed
        if seq_extract is None and self.mode == 'nt':
            seq_rec = seq_rec[::-1]  # Reverse to change phred if present
            seq_rec._seq = seq_rec._seq.complement()
            seq_extract = self._extract_seq(seq_rec)
        
        # Sequence recognized and extracted
        if seq_extract is not None:  
            self.seq = seq_extract._seq
            self.phred = seq_extract._per_letter_annotations.get('phred_quality', [])         
        else:
            self.seq = Seq('-')
            self.phred = []

        self.rfs = self._extract_rfs()
        self.quality = self._get_quality()
        # self.iso_point = self._get_pI()

    @classmethod
    def seq(cls, seq: str):
        return cls(SeqRecord(seq))

    def _select_mode(self, seq_rec: SeqRecord) -> str:
        seq_set = set(seq_rec)
        nucleotides = set('ATGCU')
        aminoacids = set('ARNDCQEGHILKMNFPOSUTWYVBZXJ*')

        if seq_set.issubset(nucleotides):
            return 'nt'
        elif seq_set.issubset(aminoacids):
            return 'aa'
        else:
            raise ValueError()

    def _load_rfs(self) -> dict:

        def _translate_to_aa(rf_dict):
            for rf, params in rf_dict.items():
                for key, value in params.items():
                    if isinstance(value, str):
                        rf_dict[rf][key] = str(Seq(value).translate(table=supE_codon_table))
                    elif isinstance(value, int):
                        rf_dict[rf][key] = value // 3
            return rf_dict
        
        file_dir = path.dirname(path.abspath(__file__))
        file_name = 'randomized_fragments.json'
        file_path = path.join(file_dir, file_name)
        
        with open(file_path, 'r') as file:
            rf_dicts = load(file)

        rf_dict = rf_dicts.get(self.lib)

        if self.mode == 'aa':
            return _translate_to_aa(rf_dict)
        return rf_dict

    def _extract_seq(self, seq_rec: SeqRecord) -> SeqRecord:

        # Remove non-letter characters except * (STOP)
        # def seq_cleanup(seq):
        #     chars = ascii_letters + '*'
        #     return ''.join([x for x in seq if x in chars])

        # seq = seq_cleanup(seq)

        # TODO Extend tag lists/objects and create algorithm
        pelB = Seq('GCGGCCCAGCCGGCCATGGCG')
        his_tag = Seq('GGCCCGGGAGGCCAACACCATCACCACCATCAT')
        myc_tag = Seq('GAACAAAAACTCATCTCAGAAGAGGATCTG')

        if self.mode == 'aa':
            pelB = pelB.translate(table=supE_codon_table)
            his_tag = his_tag.translate(table=supE_codon_table)
            myc_tag = myc_tag.translate(table=supE_codon_table)

        # Get positions of start and end regions
        pelB_i = self._get_positions(pelB, seq_rec._seq)[1]
        his_tag_i = self._get_positions(his_tag, seq_rec._seq)[0]
        myc_tag_i = self._get_positions(myc_tag, seq_rec._seq)[0]

        end_i = max(his_tag_i, myc_tag_i)

        # Add length if sequence found
        if end_i > -1:
            end_i += len(his_tag)  # TODO myc_tag

        # Extract binder sequence
        seq_extract = seq_rec[pelB_i:end_i]

        # If any sequence extracted
        if seq_extract._seq:
            return seq_extract
        else:
            return None

    def _get_positions(self, target, seq=None) -> tuple:
        
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

        alignments = aligner.align(seq, target)
        target_alignment = str(alignments[0]).split()[1]

        try:
            alignments = aligner.align(str(seq), target)
            target_alignment = str(alignments[0]).split()
            aligned_i = alignments.alignment.aligned[0][0]

            mm_condition = target_alignment.count('.') <= self.mm_lmt
            len_condition = len(range(*aligned_i)) == len(target)

            # print(target_alignment)
            # print(mm_condition, len_condition)

            if len_condition and mm_condition:
                return aligned_i
            return (-1, -1)
        except AttributeError:
            return (-1, -1)
        except IndexError:
            return (-1, -1)


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
                    extract_seq = extract_seq.translate(table=supE_codon_table)

                if pos[0] > last_pos[1]:
                    rf_extracted.append(extract_seq)
                    last_pos = pos
                else:
                    rf_extracted.append(Seq('-'))

            return rf_extracted
            
        else:
            return [Seq('-')] * len(self.rf_dict.values())


    def _lib_recognition(self) -> str:

        lib_score = {lib: 0 for lib in lib_frames.keys()}

        for lib, frames in lib_frames.items():
            for frame in frames:
                if self.mode == 'aa':
                    frame = str(Seq(frame).translate(table=supE_codon_table))
                if frame in self.seq:
                    lib_score[lib] += 1

        lib_sorted = sorted(list(lib_score.keys()),
                            key=lambda x: lib_score[x],
                            reverse=True)
        top_lib, second_lib, *_ = lib_sorted

        if lib_score[top_lib] > lib_score[second_lib]:
            return top_lib
        return 'Not recognized'

    def _get_quality(self) -> tuple:
        if self.phred:
            phred_hq = round(sum(x > self.hq_lmt for x in self.phred) / len(self.phred), 3)
            phred_min = min(self.phred)
            phred_min_i = self.phred.index(phred_min) + 1
            return (phred_hq, phred_min, f'({phred_min_i})')
        return ()

    def analyze_binder(self) -> tuple[Result, tuple]:

        def _detect_mutations() -> list:
            if self.lib != 'Not recognized':
                mutations = []
                frames = lib_frames[self.lib]
                
                if self.mode == 'aa':
                    frames = list(map(lambda x: str(Seq(x).translate(table=supE_codon_table)), frames))
                
                for i, frame in enumerate(frames, 1):
                    mutations.append(i if frame not in self.seq else '')

                return mutations
            return ['-'] * 8


        def _count_stop_codons() -> list:
            if self.mode == 'nt':
                codons = Counter([str(self.seq[i:i+3])
                                for i in range(0, len(self.seq), 3)])
                amber = codons['TAG']
                ochre = codons['TAA']
                opal = codons['TGA']
                return [n if n else '' for n in [amber, ochre, opal]]
            return ['-', '-', '-']

        # Sequence recognized
        if self.seq is not None:

            # Library detection and mutation search
            frame_muts = _detect_mutations()

            # TODO move to class
            if self.mode == 'nt':
                seq_aa = str(self.seq.translate(table=supE_codon_table))
            else:
                seq_aa = str(self.seq)

            # Stop codon search
            stop_muts = _count_stop_codons()

            # Gather results
            result = Result(self.id,
                            self.rfs,
                            self.lib,
                            self.quality,
                            frame_muts,
                            stop_muts,
                            seq_aa)

        # No sequence recognized
        else:
            result = Result(
                self.id,  # Sequence ID
                [Seq('-')] * 6,   # RFs
                'Not recognized',  # Library
                ('-', '-', '-'),  # Phred HQ, score, and position
                ['-'] * 8,   # Frame mutations
                [''] * 3,  # Stop mutations
                '-',  # Aminoacid sequence
            )
        return result


    def get_pI(self):

        # Calculates charge of aa at given pH based on pKa
        def charge(ph, pka):
            return 1 / (10 ** (ph-pka)+1)

        rfs_merged = ''.join(map(str, self.rfs)).upper()

        aa_count = Counter(rfs_merged)

        # Set range of pI search and start pH
        pH_min, pH, pH_max = 2, 8, 14.0

        base_pka = {'K': 10.00,
                    'R': 12.00,
                    'H':  5.98}
        acid_pka = {'D':  4.05,
                    'E':  4.45,
                    'C':  9.00,
                    'Y': 10.00}
        n_term_dict = {'A': 7.59,
                       'M': 7.00,
                       'S': 6.93,
                       'P': 8.36,
                       'T': 6.82,
                       'V': 7.44,
                       'E': 7.70}
        c_term_dict = {'D': 4.55,
                       'E': 4.75}

        n_pka = n_term_dict.get(rfs_merged[0], 7.7)
        c_pka = c_term_dict.get(rfs_merged[-1], 3.55)

        p_charge = sum([charge(pH, base_pka[aa]) * aa_count[aa] for aa in base_pka.keys()]) + charge(pH, n_pka)
        n_charge = sum([charge(acid_pka[aa], pH) * aa_count[aa] for aa in acid_pka.keys()]) + charge(c_pka, pH)
        net_charge = p_charge - n_charge

        while abs(net_charge) > 0.001:
            if net_charge > 0:
                pH_min = pH
            else:
                pH_max = pH

            pH = (pH_min + pH_max) / 2

            p_charge = sum([charge(pH, base_pka[aa]) * aa_count[aa] for aa in base_pka.keys()]) + charge(pH, n_pka)
            n_charge = sum([charge(acid_pka[aa], pH) * aa_count[aa] for aa in acid_pka.keys()]) + charge(c_pka, pH)
            net_charge = p_charge - n_charge

        return round(pH, 2)

    def check_unique(self, file) -> bool:
        # unique = True/False
        # return unique
        raise NotImplementedError


class RF:
    def __init__(self, params: dict, mode) -> None:
    
        self.mode = mode
        self.start = params['start']
        self.st_dis = params['start_dis']
        self.end = params['end']
        self.end_dis = params['end_dis']
        self.lmt = params['limit']
        self.end_pos = -1
        if mode == 'nt':
            self.mm_lmt = 7
        elif mode == 'aa':
            self.mm_lmt = 3

    def _get_positions(self, target, seq) -> tuple:
        aligner = Align.PairwiseAligner()
        aligner.query_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_open_gap_score = -100
        aligner.target_internal_extend_gap_score = -100
        aligner.mismatch_score = -0.5

        alignments = aligner.align(str(seq), target)
        target_alignment = str(alignments[0]).split()

        try:
            alignments = aligner.align(str(seq), target)
            target_alignment = str(alignments[0]).split()[1]
            aligned_i = alignments.alignment.aligned[0][0]

            mm_condition = target_alignment.count('.') <= self.mm_lmt
            len_condition = len(range(*aligned_i)) == len(target)

            if len_condition and mm_condition:
                return aligned_i
            return (-1, -1)
        
        except AttributeError:
            return (-1, -1)
        except IndexError:
            return (-1, -1)

    # Extract RF sequence
    def extract(self, seq) -> tuple:
        
        def _check_positions() -> bool:
            conditions = [
                self.st_pos >= 0,
                self.end_pos > self.st_pos,
                self.st_pos + self.st_dis < self.end_pos - self.end_dis,
                self.lmt >= len(seq[self.rf_start:self.rf_end])
            ]
            # print(conditions)
            return all(conditions)
        
        self.st_pos = self._get_positions(self.start, seq)[1]
        self.end_pos = self._get_positions(self.end, seq)[0]
        
        self.rf_start = self.st_pos + self.st_dis
        self.rf_end = self.end_pos - self.end_dis

        # All positions found -> adjust distance
        if _check_positions():
            seq = Seq(seq[self.rf_start:self.rf_end])
            return (seq, self.rf_start, self.rf_end)

        else:
            return (Seq('-'), 0, 0)


class SeqLiab:
    
    def __init__(self, rfs: list) -> None:
        
        if isinstance(rfs, str):
            rfs = rfs.split('_')

        self.sl_score = 0
        self.sl_count = defaultdict(int)
        self.sl_matches = defaultdict(list)

        for sl in seq_liab_db:
            sl_pattern = re.compile(sl['regex'])
            for i, rf in enumerate(rfs, 1):
                sl_matches = sl_pattern.findall(rf)
                if sl_matches:
                    self.sl_score += sl['score'] * len(sl_matches)
                    self.sl_count[sl['name']] += len(sl_matches)
                    for match in sl_matches:
                        self.sl_matches[sl['name']] += [(i, match)]

        # TODO return (i, sl_match) for every match

    def get_score(self):
        return self.sl_score

    def get_summary(self):
        return [(liab, self.sl_count[liab], self.sl_matches[liab]) for liab in self.sl_count]


if __name__ == '__main__':
    pass