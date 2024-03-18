# from DNAToolkit import *
from biobasic import *

test = biobasic()
print(test.random_nucs(35, 'DNA'))
print(test.seq_info())


# print(translation(rndSEQ, 2))

# print(codon_usage_ind(rndSEQ, 'T'))
# print(codon_usage_total(rndSEQ))

# #print(six_pack(rndSEQ))

# orf = translation(rndSEQ, 2)
# print(single_cds_find(orf))

# print(cds_sixpack_find(rndSEQ))

# print(cds_sixpack_scan(rndSEQ, 0, 120))