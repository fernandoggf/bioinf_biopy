#---------------------------------------------[DOC]------------------------------------------------------------------#
#   Methods to handle basic bioinformatic tasks with string in classes.
#   The main idea is from rebelScience.
#   Created by: Fernando G. G. Figueroa
#   Date: 17/03/24
#
#---------------------------------------------------------------------------------------------------------------------#
#----------------------------------------[-libs & metods-]------------------------------------------------------------#
import random
from bioglobal import *
#-------------------------------------------[-classes-]---------------------------------------------------------------#
class biobasic:
    """DNA strings handling class, default ATCG and DNA as type."""
    def __init__(self, seq: str ='ATCG', seq_type: str ='DNA', seq_label: str ='No Label'):
        self.seq = seq.upper()
        self.seq_type = seq_type
        self.seq_label = seq_label
        self.is_valid = self.__validate_seq()
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
        