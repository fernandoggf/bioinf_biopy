from DNAToolkit import *

# Toolkit Validations
rand = 'ACCTGSEEF'
print(validate_sequence(rand))
rand = 'TtgAC'
print(validate_sequence(rand))
print(random_nucs(33))
print(nuc_frequency(random_nucs(18)))
print(nuc_freq_coll(random_nucs(24)))
print(nuc_frequency('GTaadcFF'))
print(transcription_seq('GTACCAGT'))

# main
rndSEQ = validate_sequence(random_nucs(100))
print(reverse_complement(rndSEQ))
pretty_helix(rndSEQ)
print(gc_content(rndSEQ))
print(gc_subcontent(rndSEQ, 10)) 
