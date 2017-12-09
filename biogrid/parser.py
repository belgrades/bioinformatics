import requests
from bs4 import BeautifulSoup

url = "https://thebiogrid.org/106623/protein/homo-sapiens/adcy3.html"
r = requests.get(url)

soup = BeautifulSoup(r.text, 'html.parser')
for protein in soup.find_all('div', {'class': 'association sameOrganism'}):
    print(protein)
