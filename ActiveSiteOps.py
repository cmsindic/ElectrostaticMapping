import os
from fileProcessing import parse
import requests
from requests_html import HTMLSession


def get_uniprot_urls(url):
    ''' Locate uniprot urls from rcsb page for pdb file
    denoted by "url".
    '''
    html = requests.get(url).text
    for block in list(html.split('<')):
        if '.org/uniprot/' in block:
            uniprot_url = block.split("\"")[1]
            yield uniprot_url


def get_active_site_html(pdb):
    ''' Retrieve html of uniprot sites containing info about
    the enzyme of interest.
    '''
    generic_rcsb_url = "https://www.rcsb.org/structure/"
    generic_uniprot_url = "https://www.uniprot.org/uniprot/"
    rcsb_url = generic_rcsb_url + pdb
    headers = {
        'User-agent':
        'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/115.0'
    }
    for url in get_uniprot_urls(rcsb_url):
        session = HTMLSession()
        r = session.get(url)
        r.html.render()
        print(r.content)
