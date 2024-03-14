import requests
import xmltodict

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=3269'

base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
action = 'esumary.fcgi?'
db = 'db=gene'
term = '&id=3269'

response = requests.get(url)
if response.status_code == 200:
    data = xmltodict.parse(response.text)
    print(data)
else:
    print(f"Failed to fetch data. Status code: {response.status_code}")


