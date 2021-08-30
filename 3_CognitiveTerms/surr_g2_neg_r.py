# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 23:01:26 2021

@author: mingr
"""
from neurosynth.base.dataset import Dataset
from neurosynth import decode
import multiprocessing
import os
os.chdir(r'D:\Data\DIDA-MDD\gradient_analysis\analysis2\NeuroSynth')
dataset = Dataset.load('dataset.pkl')
terms = []
with open("g2_z_neg.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
decoder = decode.Decoder(dataset,features=terms)
def do_something(x):       
    id = '0000'+str(x) 
    os.chdir(r'D:\Data\DIDA-MDD\gradient_analysis\analysis2\NeuroSynth\surr_g2_neg')
    niiname = 'rsurr_g2_z_neg_'+id[-5:]+'.nii'
    savename = 'rsurr_g2_z_neg_'+id[-5:]+'.txt'
    decoder.decode(niiname, savename)

if __name__ == '__main__':    
    items = [x for x in range(1,10001)]
    p = multiprocessing.Pool(32)    
    p.map(do_something, items)
    