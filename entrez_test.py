import re

# API handling 
from easy_entrez import EntrezAPI
# Parsing information
from easy_entrez.parsing import xml_to_string

# Init class to use the API
entrez_api = EntrezAPI(
    'ret_toolkit',
    'fernando.ggfigueroa@icloud.com',
    return_type='json'
)

# Fetching genes for a variant from dbSNP
ret_snp = entrez_api.fetch(['rs6311'], max_results=1, database='snp')
print(ret_snp)           # <EntrezResponse status=200 for FetchQuery ['rs6311'] in snp>
print(ret_snp.data[0])   # <Element '{https://www.ncbi.nlm.nih.gov/SNP/docsum}DocumentSummary' at 0x108d22930>

ret_snp = xml_to_string(ret_snp.data[0])

# TODO: Find the genes in the xml
genes = [name for name in re.findall('<ns0:GENE_E>', ret_snp)]
print(genes)
