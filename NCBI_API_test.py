
import os
import sys 
import pandas as pd
import numpy as np
import requests
import xmltodict
from xml.etree import ElementTree
import csv

# main doc https://www.ncbi.nlm.nih.gov/books/NBK25500/


#project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
path = r"C:\Users\FernandoGonzálezGarc\Documents\bioinf\bioinf_py\bioinf_biopy\entrez_databases.tsv"
ncbi_dbs = pd.read_table(path)
entrez_dbs = ncbi_dbs['E-utility Database Name'].tolist()


actions = {'search':'esearch.fcgi?',
           'summary':'esummary.fcgi?',
           'full_download':'efetch.fcgi?',
           'linking':'elink.fcgi?',
           'global':'egquery.fcgi?'}

''' Syntax of actions
search: esearch.fcgi?db=<database>&term=<query>
             - Input: Entrez database (&db); Any Entrez text query (&term)
             - Example:  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=science[journal]+AND+breast+cancer+AND+2008[pdat]
summary: esummary.fcgi?db=<database>&id=<uid_list>
             - Input: List of UIDs (&id); Entrez database (&db)
             - Example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=6678417,9507199,28558982,28558984,28558988,28558990
full_download: efetch.fcgi?db=<database>&id=<uid_list>&rettype=<retrieval_type>&retmode=<retrieval_mode>
             - Input: List of UIDs (&id); Entrez database (&db); Retrieval type (&rettype); Retrieval mode (&retmode)
             - Example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=34577062,24475906&rettype=fasta&retmode=text
linking: elink.fcgi?dbfrom=<source_db>&db=<destination_db>&id=<uid_list>
             - Input: List of UIDs (&id); Source Entrez database (&dbfrom); Destination Entrez database (&db)
             - Example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=gene&id=34577062,24475906
global: egquery.fcgi?term=<query>
             - Input: Entrez text query (&term)
             - Example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi?term=mouse[orgn]

'''

class entrez_API:
    def __init__(self,
                toolname: str,
                email: str,
                ret_type: str,
                api_key: str = None,
                timeout: float = 10,
                url_base: str = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
                ):
        self.url_base = url_base
        self.toolname = toolname
        self.email = email
        self.ret_type: str
        self.api_key = api_key
        self.timeout = timeout
    
    def _args(self):
        return {
            'toolname': self.toolname,
            'email': self.email,
            'api_key': self.api_key
        }
    
    def _req(self, action: str, ncbi_database: str, id: list[str], rettype: str):
        # TODO:  validate args in
        params = self._args()

        construct = f'{actions[action]}db={ncbi_database}&id={id}&rettype={rettype}'
        print(construct)

my_api = entrez_API('tool', 'fer.com', 'fasta')

my_api._req('summary', 'genome', '1244', 'fasta')


# base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
# action = 'esumary.fcgi?'
# db = 'db=gene'
# query = '&id=3269'
# #TODO: construct the url complete

# response = requests.get(url)
# if response.status_code == 200:
#     data = xmltodict.parse(response.text)
#     print(data)
# else:
#     print(f"Failed to fetch data. Status code: {response.status_code}")


