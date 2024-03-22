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
        # TODO: funcion similar pero en lugar de espaciado, hacer por secciones o tramos p ej dividir en 3 la cadena 
        temp_res = []
        for i in range(0, len(self.seq)-k_spaces + 1, k_spaces):
            subseq = self.seq[i:i+k_spaces]
            temp_res.append(round((subseq.count('G')+subseq.count('C'))/len(subseq)*100, 2))
        return temp_res
    
    def translation(self, ORF: int):
        temp_trans = [codon_table[self.seq[pos:pos+3]] for pos in range(ORF, len(self.seq)-2, 3)]
        return temp_trans
