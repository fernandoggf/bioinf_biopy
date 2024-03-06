from DNAToolkit import *

# Validate nucleotides
rand = 'ACCTGSEEF'
print(validate_sequence(rand))
rand = 'TtgAC'
print(validate_sequence(rand))

print(random_nucs(33))

print(nuc_frequency(random_nucs(18)))

print(nuc_freq_coll(random_nucs(24)))

print(nuc_frequency('GTaadcFF'))