#---------------------------------------------[DOC]------------------------------------------------------------------#
#   This is a library that works for bioinformatics tasks, the main idea is from rebelScience, and adjusted for 
#   personal proposes. This toolkit has been modified through branchs and constant code. Each function has its own
#   description.
#   Created by: Fernando G. G. Figueroa
#   Date: 05/03/24
#   
#
#---------------------------------------------------------------------------------------------------------------------#
#----------------------------------------[-libs & metods-]------------------------------------------------------------#
import random
import collections 
from collections import Counter
from collections import defaultdict
#-------------------------------------------[-globals-]---------------------------------------------------------------#
nucleotides = ['A', 'C', 'G', 'T']
reverse_nucleotides = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
codon_table = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}
#-------------------------------------------[-precode-]---------------------------------------------------------------#
amino_table = defaultdict(list)
for codon, aa in codon_table.items():
    amino_table[aa].append(codon)
#------------------------------------------[-functions-]--------------------------------------------------------------#
## TODO: pretty printing, colors in console printing (video part 2), documenting every function

def validate_sequence(seq: str):
    temp_seq = seq.upper()
    for nucs in temp_seq:
        if nucs not in nucleotides:
            return False
    return temp_seq

def seq_len(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        return len(temp_seq)
    else:
        return False

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
    temp_seq = validate_sequence(seq)
    if temp_seq:
        return dict(collections.Counter(seq))
    else:
        return False

def transcription_seq(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        return seq.replace('T', 'U')
    else:
        return False
    
def reverse_complement(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_reverse = ([reverse_nucleotides[nuc] for nuc in seq])
        temp_reverse = ''.join(temp_reverse[::1])
        return temp_reverse
    else:
        return False
    
def pretty_helix(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_reverse = reverse_complement(temp_seq)
        print(f"5' {temp_seq} 3'")
        print(f"   {''.join(['|' for char in range(len(temp_seq))])}")
        print(f"3' {temp_reverse} 5'")
    else:
        return False
    
def gc_content(seq: int):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        return round((temp_seq.count('G')+temp_seq.count('C'))/len(temp_seq)*100, 3)
    else:
        return False
    
def gc_subcontent(seq: str, k_spaces: int):
    temp_seq = validate_sequence(seq)
    temp_result = []
    if temp_seq:
        for i in range(0, len(temp_seq)-k_spaces + 1, k_spaces):
            subseq = seq[i:i+k_spaces]
            temp_result.append(gc_content(subseq))
        return temp_result
    else:
        return False

def translation(seq: str, ORF: int):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_trans = [codon_table[temp_seq[pos:pos+3]] for pos in range(ORF, len(temp_seq)-2, 3)]
        return temp_trans
    else:
        return False
    
def codon_usage_ind(seq: str, aminoacid: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_list = []
        for i in range(0, len(temp_seq)-2, 3):
            if codon_table[temp_seq[i:i+3]] == aminoacid:
                temp_list.append(temp_seq[i:i+3])
        temp_dict = dict(Counter((temp_list)))
        total_weight = sum(temp_dict.values())
        for seqs in temp_dict:
            temp_dict[seqs] = round(temp_dict[seqs] / total_weight, 2)
        return temp_dict
    else:
        return False

def codon_usage_total(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_list = []
        for i in range(0, len(temp_seq)-2, 3):
            if codon_table[temp_seq[i:i+3]] in amino_table:
                temp_list.append(temp_seq[i:i+3])
        temp_dict = dict(Counter((temp_list)))
        total_weight = sum(temp_dict.values())
        for seqs in temp_dict:
            temp_dict[seqs] = round(temp_dict[seqs] / total_weight, 3)
        return temp_dict
    else:
        return False
    
def six_pack(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        frames = []
        for i in range(0,3):
            frames.append(translation(temp_seq, i))
        for i in range(0,3):
            frames.append(translation(reverse_complement(temp_seq), i))
        return frames
    else:
        return False
 
def single_cds_find(aa_seq: str):
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

def cds_sixpack_find(seq: str):
    temp_seq = validate_sequence(seq)
    if temp_seq:
        temp_frames = six_pack(temp_seq)
        six_pack_cds = []
        for orf in temp_frames:
            six_pack_cds.append(single_cds_find(orf))
        return six_pack_cds
    else:
        return False
    
def cds_sixpack_scan(seq: str, start_pos: int, end_pos: int):
    if end_pos > start_pos:
        temp_seq = validate_sequence(seq[start_pos:end_pos])
        if temp_seq:
            return cds_sixpack_find(temp_seq)
        else:
            return False
    else:
        return False

# Next ncbi python api, make class functions