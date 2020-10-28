# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 21:08:39 2020

@author: jinqi
"""

import numpy as np

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
        
            