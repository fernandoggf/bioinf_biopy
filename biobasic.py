#---------------------------------------------[DOC]------------------------------------------------------------------#
#   Methods to handle basic bioinformatic tasks with string in classes.
#   The main idea is from rebelScience.
#   Created by: Fernando G. G. Figueroa
#   Date: 17/03/24
#
#---------------------------------------------------------------------------------------------------------------------#
#----------------------------------------[-libs & metods-]------------------------------------------------------------#
import random
from collections import Counter
from bioglobal import *
#-------------------------------------------[-classes-]---------------------------------------------------------------#
# TODO: add docstring to every function
class biobasic:
    """DNA strings handling class, default ATCG and DNA as type."""
    def __init__(self, seq: str ='ATCG', seq_type: str ='DNA', seq_label: str ='No Label'):
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.seq_label = seq_label
        self.is_valid = self.__validate_seq()
        self.seq_rna = None
        self.second_strand = None
        assert self.is_valid, f'Provided data: {seq}, seems not be correct in the context of this method.'

    def __validate_seq(self):
        return set(nucleotides).issuperset(self.seq)
    
    def seq_info(self):
        return f'[Type]: {self.seq_type}\n[Label]: {self.seq_label}\n[Sequence]: {self.seq}\n[Length]: {len(self.seq)}'
    
    def seq_length(self):
        return len(self.seq)
    
    def random_nucs(self, len_seq: int = 10, seq_type:str ='DNA'):
        seq = ''.join([random.choice(nucleotides) for nuc in range(len_seq)])
        self.__init__(seq, seq_type, 'Randomly Generated Sequence (RGC)')
        return seq
    
    def nuc_frequency(self):
        return dict(Counter(self.seq))

    def transcription(self):
        self.seq_rna = self.seq.replace('T', 'U')
        return self.seq_rna

    def reverse_complement(self):
        temp_reverse  =  ([reverse_nucleotides[nuc] for nuc in self.seq])
        self.second_strand  = ''.join(temp_reverse[::1])
        return self.second_strand
    
    def gc_content(self):
        return round((self.seq.count('G')+self.seq.count('C'))/len(self.seq)*100, 2)
    
    def gc_subcontent(self, k_spaces: int):
        temp_res = []
        for i in range(0, len(self.seq)-k_spaces + 1, k_spaces):
            subseq = self.seq[i:i+k_spaces]
            temp_res.append(round((subseq.count('G')+subseq.count('C'))/len(subseq)*100, 2))
        return temp_res
    
    def gc_subcon_sec(self, parts: int):
        temp_res = []
        part_size = len(self.seq) // parts
        extra_chars = len(self.seq) % parts
        subseqs = []
        start_index = 0
        for i in range(parts):
            part_length = part_size + (1 if i < extra_chars else 0)
            subseqs.append(self.seq[start_index:start_index+part_length])
            start_index += part_length
        for seqs in subseqs:
            temp_res.append(round((seqs.count('G')+seqs.count('C'))/len(seqs)*100, 2))
        return temp_res

    def translation(self, ORF: int):
        temp_trans = [codon_table[self.seq[pos:pos+3]] for pos in range(ORF, len(self.seq)-2, 3)]
        return temp_trans
    
    def codon_usage(self, total: bool, aminoacid: str=None):
        temp_list = []
        for i in range(0, len(self.seq)-2, 3):
            if total:
                if codon_table[self.seq[i:i+3]] in amino_table:
                    temp_list.append(self.seq[i:i+3])
            else:
                if codon_table[self.seq[i:i+3]] == aminoacid:
                    temp_list.append(self.seq[i:i+3])
        temp_dict = dict(Counter((temp_list)))
        total_weight = sum(temp_dict.values())
        for seqs in temp_dict:
            temp_dict[seqs] = round(temp_dict[seqs] / total_weight, 2)
        return temp_dict
    
    def six_pack(self):
        frames = []
        for i in range(0,3):
            frames.append(self.translation(i))
        temp_seq = biobasic(self.reverse_complement())
        for i in range(0,3):
            frames.append(temp_seq.translation(i))
        del temp_seq
        return frames
    
    def single_cds_find(self, aa_seq: str):
        temp_cds = []
        cds = []
        for aa in aa_seq:
            if aa == '*':
                if temp_cds:
                    cds.append(temp_cds)
                    temp_cds = []
            else: 
                if aa == 'M':
                    temp_cds.append('')
                for i in range(len(temp_cds)):
                    temp_cds[i] += aa 
        return cds
    
    def cds_sixpack_find(self):
        temp_frames = self.six_pack()
        six_pack_cds = []
        for orf in temp_frames:
            six_pack_cds.append(self.single_cds_find(orf))
        return six_pack_cds
    # TODO; digest result: enlist the largest to smallest protein and
    # says from what orf it came

    def cds_sixpack_scan(self, start_pos: int, end_pos: int):
        if end_pos > start_pos:
            temp_seq = biobasic(self.seq[start_pos:end_pos])
            return temp_seq.cds_sixpack_find()
        else:
            return False
    # TODO: well document how this scans works bcs is from six on six