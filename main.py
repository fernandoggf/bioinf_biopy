from DNAToolkit import *

# main
rndSEQ = validate_sequence('ATAAGCATTCTCTGATGAAGCTTTGCGCGGTCTTGCGTACGGGAGTAGCACAGCC')
print(reverse_complement(rndSEQ))
pretty_helix(rndSEQ)
print(gc_content(rndSEQ))
print(gc_subcontent(rndSEQ, 5)) 

print(translation(rndSEQ, 0))

print(codon_usage_ind(rndSEQ, 'T'))
print(codon_usage_total(rndSEQ))

print(six_pack(rndSEQ))