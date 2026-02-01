# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 19:49:36 2023

@author: qihan
"""

import pandas as pd

def read_pickle_file(file):
    pickle_data = pd.read_pickle(file)
    return pickle_data


'''
import pickle

with open('N5236.pkl', 'rb') as f:
    data = pickle.load(f)
    
'''