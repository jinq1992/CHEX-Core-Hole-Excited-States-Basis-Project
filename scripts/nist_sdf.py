# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:44:33 2019

@author: jinqi
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Download all information such as UV spectrum, mol file available from NIST Chemistry Webbook."""
import os
import re

import requests
from bs4 import BeautifulSoup


NIST_URL = 'http://webbook.nist.gov/cgi/cbook.cgi'
EXACT_RE = re.compile('/cgi/cbook.cgi\?GetInChI=(.*?)$')
ID_RE = re.compile('/cgi/cbook.cgi\?ID=(.*?)&')
JDX_PATH = 'nist/jdx/'
MOL_PATH = 'nist/mol/'
SDF_PATH = 'nist/sdf/'
##3dsdf file https://webbook.nist.gov/cgi/cbook.cgi?Str3File=C5314830
##2dmol file https://webbook.nist.gov/cgi/cbook.cgi?Str2File=C5314830

def search_nist_inchi(inchi):
    """Search NIST using the specified InChI or InChIKey query and return the matching NIST ID."""
    print('Searching: %s' % inchi)
    response = requests.get(NIST_URL, params={'InChI': inchi, 'Units': 'SI'})
    soup = BeautifulSoup(response.text)
    idlink = soup.find('a', href=EXACT_RE)
    if idlink:
        nistid = re.match(EXACT_RE, idlink['href']).group(1)
        print('Result: %s' % nistid)
        return nistid
    # If no match, there is a list of similar species
    # if not ids:
    #     ids = [re.match(ID_RE, link['href']).group(1) for link in soup('a', href=ID_RE)]


def search_nist_formula(formula, allow_other=False, allow_extra=False, match_isotopes=False, exclude_ions=False, has_uv=False):
    """Search NIST using the specified formula query and return the matching NIST IDs."""
    print('Searching: %s' % formula)
    params = {'Formula': formula, 'Units': 'SI'}
    if allow_other:
        params['AllowOther'] = 'on'
    if allow_extra:
        params['AllowExtra'] = 'on'
    if match_isotopes:
        params['MatchIso'] = 'on'
    if exclude_ions:
        params['NoIon'] = 'on'
    if has_uv:
        params['cUV'] = 'on'
    response = requests.get(NIST_URL, params=params)
    soup = BeautifulSoup(response.text)
    ids = [re.match(ID_RE, link['href']).group(1) for link in soup('a', href=ID_RE)]
    print('Result: %s' % ids)
    return ids


#def get_jdx(nistid, stype='UVVis'):
#    """Download jdx file for the specified NIST ID, unless already downloaded."""
#    filepath = os.path.join(JDX_PATH, '%s-%s.jdx' % (nistid, stype))
#    if os.path.isfile(filepath):
#        print('%s %s: Already exists at %s' % (nistid, stype, filepath))
#        return
#    print('%s %s: Downloading' % (nistid, stype))
#    response = requests.get(NIST_URL, params={'JCAMP': nistid, 'Type': stype, 'Index': 0})
#    if response.text == '##TITLE=Spectrum not found.\n##END=\n':
#        print('%s %s: Spectrum not found' % (nistid, stype))
#        return
#    print('Saving %s' % filepath)
#    with open(filepath, 'w') as file:
#        file.write(response.content)


def get_sdf(nistid):
    """Download mol file for the specified NIST ID, unless already downloaded."""
    filepath = os.path.join(SDF_PATH, '%s.mol' % nistid)
    if os.path.isfile(filepath):
        print('%s: Already exists at %s' % (nistid, filepath))
        return
    print('%s: Downloading sdf' % nistid)
    response = requests.get(NIST_URL, params={'Str3File': nistid})
    if response.text == 'NIST    12121112142D 1   1.00000     0.00000\nCopyright by the U.S. Sec. Commerce on behalf of U.S.A. All rights reserved.\n0  0  0     0  0              1 V2000\nM  END\n':
        print('%s: SDF not found' % nistid)
        return
    print('Saving %s' % filepath)
    with open(filepath, 'wb') as file:
        file.write(response.content)


#def get_all_uvvis():
#    """Search NIST for all structures with UV-Vis spectra and download a JDX file for each."""
#    # Each search is limited to 400 results
#    # So we search by formula, allowing additional elements not specified in formula: C, CC, CCC, CCCC, etc.
#    for i in range(1, 100):
#        ids = search_nist_formula('C%s' % i, allow_other=True, exclude_ions=True, has_uv=True)
#        print('%s spectra found' % len(ids))
#        for nistid in ids:
#            get_mol(nistid)
#            get_jdx(nistid, stype='UVVis')


#if __name__ == '__main__':
    #nistid = search_nist_inchi('ZYGHJZDHTFUPRJ-UHFFFAOYSA-N')
    #get_jdx(nistid)
    #get_mol(nistid)
    #search_nist_formula('C20', allow_other=True, exclude_ions=True, has_uv=True)
 #   get_all_uvvis()
    
####the following directly address XPS_accuracy project
    
import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
#df = pd.read_csv('file:///C:/Users/jinqi/XPS_benchmark.csv')
#df['nistid'] = 'NaN'
#for index,row in df.iterrows():
#    formula = row['Formula']
#    nist_id = search_nist_formula(formula)
#    if nist_id == []:
#        df['nistid'][index] = 'NaN'
#    if nist_id != []:
#        df['nistid'][index] = nist_id[0]
#
#df.to_csv('new_XPS_benchmark.csv')
df = pd.read_csv('file:///C:/Users/jinqi/new_XPS_benchmark.csv')
for index,row in df.iterrows():
    if row['nistid'] =='NaN':
        continue
    if row['nistid'] !='NaN':
        get_sdf(row['nistid'])
        