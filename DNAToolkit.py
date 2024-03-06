# TODO: docstring of this file 
import random
import collections

# My own toolkit to work with
nucleotides = ['A', 'C', 'G', 'T']

def validate_sequence(seq: str):
    temp_seq = seq.upper()
    for nucs in temp_seq:
        if nucs not in nucleotides:
            return False
    return temp_seq

def random_nucs(len_seq: int):
    rand_seq = ''.join([random.choice(nucleotides) for nuc in range(len_seq)])
    return rand_seq

def nuc_frequency(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_freq = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
        for nucs in temp_seq:
            temp_freq[nucs] += 1
        return temp_freq
    else:
        return False

def nuc_freq_coll(seq: str):
    return dict(collections.Counter(seq))