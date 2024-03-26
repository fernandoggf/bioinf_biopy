# from DNAToolkit import *
from biobasic import *

test = biobasic()
print(test.random_nucs(200, 'DNA'))
print(test.seq_info())

print(test.translation(0))
print(test.gc_subcon_sec(5))

print(test.codon_usage(total=True))

for frames in test.six_pack():
    print(frames)

print(test.cds_sixpack_find())

print(test.cds_sixpack_scan(20, 80))
