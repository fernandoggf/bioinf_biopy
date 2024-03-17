
import os
import sys 
import csv
import time
import pandas as pd
import numpy as np
import requests
import xmltodict
from xml.etree import ElementTree
from time import sleep


# main doc https://www.ncbi.nlm.nih.gov/books/NBK25500/

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
ncbi_dbs = pd.read_table(os.path.join(project_dir, 'bioinf_biopy', 'entrez_databases.tsv'))
entrez_dbs = ncbi_dbs['E-utility Database Name'].tolist()

actions = {'search':'esearch.fcgi?',
           'summary':'esummary.fcgi?',
           'full_download':'efetch.fcgi?',
           'linking':'elink.fcgi?',
           'global':'egquery.fcgi?'}

retrieve_type = ['fasta', 'gbwithparts', 'gb', 'uilist', 'abstract', 'pdb']

retrieve_mode = ['xml', 'json', 'html']

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
                #ret_type: str,
                api_key: str = None,
                timeout: float = 10,
                url_base: str = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
                min_time_gdl: float = 0.338
                ):
        self.url_base = url_base
        self.toolname = toolname
        self.email = email
        self.ret_type: str
        self.api_key = api_key
        self.timeout = timeout
        self.last_req = None
    
    def _args(self):
        return {
            'toolname': self.toolname,
            'email': self.email,
            'api_key': self.api_key
        }
    
    # TODO: input parameter termo to access like this example https://www.ncbi.nlm.nih.gov/nuccore?term=homo%20sapiens%5BOrganism%5D
    def _req(self, action: str, ncbi_database: str, id: list[str], rettype: str, retmode: str = 'xml'):
        param = self._args()
        # Validate inputs
        if action not in actions:
            print('User error EU1: invalid action to do')
            sys.exit()
        if ncbi_database not in entrez_dbs:
            print('User error EU2: invalid database to access')
            sys.exit() 
        if rettype not in retrieve_type:
            print('User error EU3: invalid retrieve type to pass')
            sys.exit()
        if retmode not in retrieve_mode:
            print('User error EU4: invalid retrieve mode to receive')
            sys.exit()
        
        # Construction of url
        construct = f'{self.url_base}{actions[action]}db={ncbi_database}&id={id}&rettype={rettype}&retmode={retmode}'
        print(construct)

        # time validation to handle guidelines
        current_time = time.time()
        if self.last_req != None:
            time_btwn = current_time - self.last_req
            if time_btwn < self.min_time_gdl:
                time_to_wait = self.min_time_gdl - time_btwn
                sleep(time_to_wait)
        self.last_req = current_time

        # request handle
        try:
            response = requests.get(construct, params=param, timeout=self.timeout)
            if response.status_code == 200:
                # TODO: make a list the reponse
                print(response.text)
            else:
                print(f"Request failed with status code {response.status_code}")
        except requests.exceptions.Timeout:
            print("Timeout error: The request took too long to complete.")
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")


def main():
    my_api = entrez_API('ret_tool', 'fernando.ggfigueroa@icloud.com')

    my_api._req('summary', 'nuccore', '2244', 'fasta', 'xml')

if __name__ == "__main__":
    main()


