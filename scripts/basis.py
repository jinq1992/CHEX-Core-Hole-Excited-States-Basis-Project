# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:55:08 2020

@author: jinq
"""
import sys
print(sys.path)

'''this code will read basis sets established in the literature and plot the electronic wave function'''
import matplotlib.pyplot as plt
import numpy as np
import requests
from bs4 import BeautifulSoup
import urllib.request
import urllib.parse
import re
import pandas as pd

print ('This code will read basis sets established in the literature and plot the electronic wave function. \
Functions include: BSE(basis_set_choice,element_choice), plot_qchemto1s(content1s), plot_qchemto2s(content2s), plot_qchemto2p(content2p), scale(content,scaling_factor) ect.  ')

C1s1 = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C1s1.csv', delimiter=',',names = ['rC1s1','rRC1s1'])
rC_1s_gs = C1s1['rC1s1'].to_numpy()
x = rC_1s_gs

#The code will first print out the basis set information
#Then will show which part is used to construct atomic orbital
#Then reads USPP generated AO
#Then compares constructed atomic orbital vs USPP orbital 


#basic settings and global variable
# =============================================================================
# basis_set_choice = input('suggested basis set choices include: 6-31++G**, cc-pVTZ, def2-TZVP, IGLO-II, Sapporo-TZP... what is your basis set choice? ')
# print('you have chosen: ' + basis_set_choice)
# element_choice = input('suggested choice of element include: C, N, H ect. for this function only single element is allowed ... what is your choice? ')
# print('you have chosen: ' + element_choice)
# =============================================================================
##end of user input basic seetings and global variable

def BSE(basis_set_choice,element_choice):
    '''this function will get the basis set on BSE, and return the qchem basis format\
    source is from https://www.basissetexchange.org'''
    # now get the basis set result for q-chem type input from BS Exchange site
    url =  'https://www.basissetexchange.org/'
    #example https://www.basissetexchange.org/basis/def2-tzvp/format/qchem/?version=1&elements=6
    #example https://www.basissetexchange.org/basis/cc-pvtz/format/qchem/?version=1&elements=6
    #example https://www.basissetexchange.org/basis/6-31%2B%2Bg_st__st_/format/qchem/?version=1&elements=6
    #example https://www.basissetexchange.org/basis/sapporo-tzp/format/qchem/?version=1&elements=6
    basis_dic = {'6-31++G**':'basis/6-31%2B%2Bg_st__st_/format/qchem/?version=1',
                 '6-31++G*':'basis/6-31%2B%2Bg_st_/format/qchem/?version=1',
                 '6-31++G':'basis/6-31%2B%2Bg/format/qchem/?version=1',
                 '6-31+G':'basis/6-31%2Bg/format/qchem/?version=1',
                 '6-31+G*':'basis/6-31%2Bg_st_/format/qchem/?version=1',
                 '3-21G':'basis/3-21g/format/qchem/?version=1',
                 'cc-PVTZ':'basis/cc-pvtz/format/qchem/?version=1',
                 'aug-cc-PVTZ':'basis/aug-cc-pvtz/format/qchem/?version=1',
                 'def2-TZVP':'basis/def2-tzvp/format/qchem/?version=1',
                 'Sapporo-TZP':'basis/sapporo-tzp/format/qchem/?version=1'}
    #only the first 3 rows of periodic table is included in the current code
    ele_dic   = {'H':'&elements=1','He':'&elements=2','Li':'&elements=3',
                 'Be':'&elements=4','B':'&elements=5','C':'&elements=6',
                 'N':'&elements=7','O':'&elements=8','F':'&elements=9',
                 'Ne':'&elements=10','Na':'&elements=11','Mg':'&elements=12',
                 'Al':'&elements=13','Si':'&elements=14','P':'&elements=15',
                 'S':'&elements=16','Cl':'&elements=17','Ar':'&elements=18'}
    new_url = url + basis_dic[basis_set_choice] + ele_dic[element_choice]
    print()
    print('you can find the basis inputs at the following website:')
    print(new_url)
    web_result = requests.get(new_url)
    if web_result.status_code != 200:
        #for other potential status code, see https://en.wikipedia.org/wiki/List_of_HTTP_status_codes
        print()
        print('error in accessing the BSE website, check again')
    if web_result.status_code == 200:
        print()
        print('successfully access to the BSE website for ' + basis_set_choice + element_choice)
        ##lastly, parse the basis set info for future usage using BeautifulSoup
        src = web_result.content #store the content as source as src
        soup = BeautifulSoup(src,'lxml')
        content = soup.find(id = 'bstext') #only getting the content
        print()
        print(soup.find(id = 'bstext')) #print the basis set XXX for element XXX
        return content
  
# =============================================================================
#run the function, as an example we run 6-31++G** for C
#BSE(basis_set_choice,element_choice)
#basis_set_choice = 'cc-PVTZ'
#element_choice   = 'C'
#bs_content = str(BSE(basis_set_choice,element_choice))
#lst = bs_content.split('\n')
##lst is the most important attribute, which is the whole basis set content
# =============================================================================

# =============================================================================
#################=======================1s orbital=============================
# =============================================================================
def counting_1s(lst):
    ''' this function will read the basis set text, count line numbers for the start and end of 1s orbital'''
    counter = 0
    for i in lst:
        counter += 1
        if 'S   ' in i:
            no_gaussian = int(i.split()[1])
            start_1s = counter
            return (start_1s, no_gaussian)           
def qchemto1s(line_start_1s,line_end_1s,lst):
    '''this function will read line numbers, and return 1s orbital content1s'''
    #processing the 1s part of the results in BSE()
    #disecting the BSE into 1s part
    print ('running function qchem21s()')
    print ('you are using the content for creating 1s xyz orbital')
    content1s = lst[line_start_1s-1:line_end_1s]
    for i in content1s:
        print (i)
    return content1s
def plot_qchemto1s(x,content1s):
    '''this function will use the content1s in qchemto1s() and plot it'''
    line_1s = content1s[1:]
    GTF = 0 #initiation of the GTO from GTF
    for i in line_1s:
        exponential = i.split()[0]
        exp_para = float(exponential.replace('D','E')) #replace the annoying notation
        contraction = i.split()[1]
        con_para = float(contraction.replace('D','E')) #and then change to float, precision is perserved
        #equation for constructing 1s type GTO 
        GTF += con_para*(2*exp_para/np.pi)**(3/4)*np.exp(-exp_para*x**2)
    #plt.plot(x, GTF, color = 'black')
    #plt.show()
    #plt.xlim(0,3)
    print('returning 1s orbital from ' + content1s[0])
    return GTF
# =============================================================================
##reminder to check out other basis set, not necessary finding S in the sentence is good
##execution of 1s part
#line_start_1s = counting_1s(lst)[0]
#line_end_1s = line_start_1s + counting_1s(lst)[1] 
#content1s = qchemto1s(line_start_1s,line_end_1s,lst)
#plot_qchemto1s(content1s)
# ============================================================================= 


# ============================================================================= 
#################=======================2s orbital=============================
# =============================================================================        
def counting_2s(lst,basis_set_choice):
    ''' this function will read the basis set text, find the start and end of 2s orbital'''
    if basis_set_choice == 'cc_PVTZ' or 'aug-cc-PVTZ':
        occurance = 3
    elif '6-31' or '3-21' in basis_set_choice:
        occurance = 4
    line_no_2s = [i for i, n in enumerate(lst) if 'S' in n][occurance] #3 for cc-pvtz, 4 for pople   
    start_2s = line_no_2s + 1
    line_2s = lst[line_no_2s]
    no_gaussian = int(line_2s.split()[1])
    return (start_2s, no_gaussian)
def qchemto2s(line_start_2s,line_end_2s,lst):
    '''this function will take line numbers and return the qchem basis 2s orbital'''
    #processing the 1s part of the results in BSE()
    #disecting the BSE into 1s part
    print ('running function qchemto2s()')
    print ('you are using the content for creating 2s xyz orbital')
    content2s = lst[line_start_2s-1:line_end_2s]
    for i in content2s:
        print (i)
    return content2s
def plot_qchemto2s(x, content2s):
    '''this function will use the content1s in qchemto2s() and plot it'''
    line_2s = content2s[1:]
    GTF = 0 #initiation of the GTO from GTF
    for i in line_2s:
        exponential = i.split()[0]
        exp_para = float(exponential.replace('D','E')) #replace the annoying notation
        contraction = i.split()[1]
        con_para = float(contraction.replace('D','E')) #and then change to float, precision is perserved
        #equation for constructing 2s type GTO 
        GTF += -con_para*(2*exp_para/np.pi)**(3/4)*np.exp(-exp_para*x**2)
#    plt.plot(x, GTF, label = 'gs_2s_' + basis_set_choice + '_' + element_choice, color = 'red')
#    plt.show()
#    plt.xlim(0,3)
#    plt.legend()
    return GTF
    print('plotting 2s orbital from ' + content2s[0])
# =============================================================================
##reminder to check out other basis set, not necessary finding SP in the sentence is good
##execution of 2s part
#line_start_2s = int(counting_2s(lst,basis_set_choice)[0])
#line_end_2s = line_start_2s + int(counting_2s(lst,basis_set_choice)[1])
#content2s = qchemto2s(line_start_2s,line_end_2s,lst)
#plot_qchemto2s(content2s)

# =============================================================================
def calculateRMSE(y1,y2):
    RMSE= np.sqrt(((y1-y2) ** 2).mean())
    print('RMSE is calculated to be ' + str(RMSE))

# ============================================================================= 
#################=======================2p orbital=============================
# =============================================================================        
def counting_2p(lst):
    ''' this function will read the basis set text, find the start and end of 2p orbital'''
    counter = 0
    for i in lst:
        counter += 1
        if 'P   ' in i:  #the keyword to look for is SP, be careful of this assumption for non-pople basis sets
            no_gaussian = int(i.split()[1])
            start_2p = counter
            return (start_2p, no_gaussian)
def qchemto2p(line_start_2p,line_end_2p,lst):
    '''this function will take line numbers and return the qchem basis 2p orbital'''
    #processing the 1s part of the results in BSE()
    #disecting the BSE into 1s part
    print ('running function qchemto2p()')
    print ('you are using the content for creating 2p xyz orbital')
    content2p = lst[line_start_2p-1:line_end_2p]
    for i in content2p:
        print (i)
    return content2p
def plot_qchemto2p(x, content2p):
    '''this function will use the content2p in qchemto2p() and plot it'''
    line_2p = content2p[1:]
    GTF = 0 #initiation of the GTO from GTF
    for i in line_2p:
        exponential = i.split()[0]
        exp_para = float(exponential.replace('D','E')) #replace the annoying notation
        contraction = i.split()[-1] #-1 for both pople and dunning
        con_para = float(contraction.replace('D','E')) #and then change to float, precision is perserved
        #equation for constructing 2s type GTO 
#        GTF += x*con_para*((2*exp_para/np.pi)**(3/4))*np.exp(-exp_para*x**2)
        GTF += x*con_para*((2*exp_para/np.pi)**(3/4))*((4/3*exp_para)**(1/2))*np.exp(-exp_para*x**2)          #1.3 is a magic number that gives a perfect fit, I don't know why
#        y = rRC_2p_gs/rC_2p_gs 
#        RMSE= np.sqrt(((y[0:120]-GTF ) ** 2).mean())
#        print(RMSE)   
#    plt.plot(x, GTF, label = 'gs_2p_' + basis_set_choice + '_' + element_choice, color = 'blue')
#    plt.show()
#    plt.xlim(0,3)
#    plt.legend()
    return GTF
    print('plotting 2p orbital from ' + content2p[0])
# =============================================================================
##reminder to check out other basis set, not necessary finding SP in the sentence is good
##execution of 2p part
#line_start_2p = int(counting_2p(lst)[0])
#line_end_2p = line_start_2p + int(counting_2p(lst)[1])
#content2p = qchemto2p(line_start_2p,line_end_2p,lst)
#plot_qchemto2p(content2p)
# =============================================================================

def uspp21s():
    '''this function will take USPP 1s orbital as input and plot it'''
C1s1 = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C1s1.csv', delimiter=',',names = ['rC1s1','rRC1s1'])
C1s_gs = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C1s_gs.csv', delimiter=',',names = ['rCgs','rRCgs'])
N1s1 = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\N1s1.csv', delimiter=',',names = ['rN1s1','rRN1s1'])
N1s_gs = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\N1s_gs.csv', delimiter=',',names = ['rNgs','rRNgs'])
rC_1s_gs = C1s_gs['rCgs'].to_numpy()
rC_1s1 = C1s1['rC1s1'].to_numpy()
rN_1s_gs = N1s_gs['rNgs'].to_numpy()
rN_1s1 = N1s1['rN1s1'].to_numpy()
rRC_1s_gs = C1s_gs['rRCgs'].to_numpy()*7.66030807/27.3371105
rRC_1s1 = C1s1['rRC1s1'].to_numpy()*7.66030807/27.3371105
rRN_1s_gs = N1s_gs['rRNgs'].to_numpy()*7.66030807/27.3371105
rRN_1s1 = N1s1['rRN1s1'].to_numpy()*7.66030807/27.3371105
#plt.plot(rC_1s_gs,rRC_1s_gs/rC_1s_gs,label = 'C1s_gs_AE',linestyle='dashed',color ='black', linewidth = 3)
#plt.plot(rC_1s1,rRC_1s1/rC_1s1,label = 'C1s1_AE',linestyle='dashed',color ='black', linewidth = 3)
#plt.legend()
    
def uspp22s():
    '''this function will take USPP 2s obrtial as input and plot it'''
C2s1 = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C2s1.csv', delimiter=',',names = ['rC2s1','rRC2s1'])
C2s_gs = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C2s_gs.csv', delimiter=',',names = ['rCgs','rRCgs'])
rC_2s_gs = C2s_gs['rCgs'].to_numpy()
rRC_2s_gs = C2s_gs['rRCgs'].to_numpy()*7.66030807/27.3371105
rC_2s1 = C2s1['rC2s1'].to_numpy()
rRC_2s1 = C2s1['rRC2s1'].to_numpy()*7.66030807/27.3371105
#plt.plot(rC_2s_gs,rRC_2s_gs/rC_2s_gs,label = 'C2s_gs_AE',linestyle='dashed',color ='red', linewidth = 3)
#plt.plot(rC_2s1,rRC_2s1/rC_2s1,label = 'C2s1_AE',linestyle='dashed',color ='red', linewidth = 3)

#plt.legend()

  
def uspp22p():
    '''this function will take USPP 2p orbital as input and plot it'''
C2p1 = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C2p1.csv', delimiter=',',names = ['rC2p1','rRC2p1'])
C2p_gs = pd.read_csv('C:\\Users\\jinqi\\Desktop\\USPP_basis\\C2p_gs.csv', delimiter=',',names = ['rCgs','rRCgs'])
rC_2p_gs = C2p_gs['rCgs'].to_numpy()
rRC_2p_gs = C2p_gs['rRCgs'].to_numpy()*7.66030807/27.3371105
rC_2p1 = C2p1['rC2p1'].to_numpy()
rRC_2p1 = C2p1['rRC2p1'].to_numpy()*7.66030807/27.3371105  
y = rRC_2p_gs/rC_2p_gs 
#plt.plot(rC_2p_gs,y[0:420],label = 'C2p_gs_AE',linestyle='dashed',color ='blue', linewidth = 3)
#plt.plot(rC_2p1,rRC_2p1/rC_2p1,label = 'C2p1_AE',linestyle='dashed',color ='blue', linewidth = 3)
#plt.legend()


############======================execution=======================
               
def scale(content, x):
    '''this function takes in basis function, and scale the exponential part accordingly to the scaling factor x'''
    lst = list()
    lst.append(content[0])
    for i in content[1:]:
            if 'D' in i:
                j = i.replace('D','E')
                new_exp = float(j.split()[0])*x
                old_con = float(j.split()[1])
            if 'D' not in i:
                new_exp = float(i.split()[0])*x
                old_con = float(i.split()[1])
            lst.append(str(new_exp) + ' ' + str(old_con))
    return(lst)


#fitted gaussians for excited states
#content1s1 = ['S   10   1.00',
#'1.485484D+03	 1.727672D-01',
#'1.138436D+03	-4.755598D-01',
#'5.160074D+02	-6.373887D-01',
#'7.469061D+02	7.491032D-01',
#'1.984919D+02	-7.533147D-01',
#'2.336357D+02	9.291254D-01',
#'6.371942D+01	-1.985973D+00',
#'6.385802D+01	2.187829D+00',
#'1.384818D+01	4.956949D-01',
#'3.758450D+00	4.567014D-01']
#plot_qchemto1s(content1s)
#
#content2s1 = ['S   10   1.00',
#'4.109563D+02   -5.977572D-02',
#'1.154028D+01	2.961690D+00',
#'2.651179D+02	1.634669D-01',
#'1.166457D+02	1.401256D-01',
#'1.699591D+02	-1.791583D-02',
#'1.030128D+02	-6.461064D-03',
#'1.517634D+02	-2.458007D-01',
#'7.822145D-01	1.225464D+01',
#'1.169499D+01	-3.218665D+00',
#'8.269691D-01	-1.137180D+01']
#plot_qchemto2s(content2s1) 
#
#content2p1 = ['P   5   1.00',
#'2.417898D+01	3.665224D-01',
#'2.415470D+01	-2.248861D-01',
#'3.113827D+00	1.307634D+00',
#'3.114298D+00	-6.980677D-01',
#'5.601885D-01	6.953299D-01']
#plot_qchemto2p(content2p1)
#USPP plot for excited states
    