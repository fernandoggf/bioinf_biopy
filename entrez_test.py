import re
# API handling 
from easy_entrez import EntrezAPI
# Parsing information
from easy_entrez.parsing import xml_to_string

import xml.etree.ElementTree as ET

# Init class to use the API
entrez_api = EntrezAPI(
    'ret_toolkit',
    'fernando.ggfigueroa@icloud.com',
    return_type='json'
)

# Fetching genes for a variant from dbSNP
ret_snp = entrez_api.fetch(['rs6311', 'rs662138'], max_results=2, database='snp')
print(ret_snp)           # <EntrezResponse status=200 for FetchQuery ['rs6311'] in snp>
print(ret_snp.data[0])   # <Element '{https://www.ncbi.nlm.nih.gov/SNP/docsum}DocumentSummary' at 0x108d22930>
print(xml_to_string(ret_snp.data[1]))

ret_snp_xml = ET.fromstring(xml_to_string(ret_snp.data[1]))

maf_elements = ret_snp_xml.findall('.//{https://www.ncbi.nlm.nih.gov/SNP/docsum}MAF')
for maf in maf_elements:
    study = maf.find('{https://www.ncbi.nlm.nih.gov/SNP/docsum}STUDY').text
    freq = maf.find('{https://www.ncbi.nlm.nih.gov/SNP/docsum}FREQ').text
    print(f"Study: {study}, Frequency: {freq}")

genes_element = ret_snp_xml.find('.//{https://www.ncbi.nlm.nih.gov/SNP/docsum}GENES')
if genes_element is not None:
    gene_name = genes_element.find('{https://www.ncbi.nlm.nih.gov/SNP/docsum}GENE_E/{https://www.ncbi.nlm.nih.gov/SNP/docsum}NAME').text
    gene_id = genes_element.find('{https://www.ncbi.nlm.nih.gov/SNP/docsum}GENE_E/{https://www.ncbi.nlm.nih.gov/SNP/docsum}GENE_ID').text
    
    print(f"Gene Name: {gene_name}")
    print(f"Gene ID: {gene_id}")
else:
    print("No GENES element found in the XML.")
