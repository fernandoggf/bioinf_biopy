from DNAToolkit import *

# main
rndSEQ = validate_sequence('ATATGAGCGTAAGCATTCATTCTCTGAACTTGCGCTGCGTATTGCGGCGTAATTCTGCTTGCCGTACTGCGGTCTTGCGTACCTTGCGGTCTTGTGCGGGAGTAGCACTTGCAGCCTGA')
print(reverse_complement(rndSEQ))
#pretty_helix(rndSEQ)
print(gc_content(rndSEQ))
print(gc_subcontent(rndSEQ, 5)) 

print(translation(rndSEQ, 2))

print(codon_usage_ind(rndSEQ, 'T'))
print(codon_usage_total(rndSEQ))

#print(six_pack(rndSEQ))

orf = translation(rndSEQ, 2)
print(single_cds_find(orf))

print(cds_sixpack_find(rndSEQ))

print(cds_sixpack_scan(rndSEQ, 0, 120))