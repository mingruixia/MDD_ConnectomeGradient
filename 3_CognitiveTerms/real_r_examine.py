# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 20:26:03 2021

@author: mingr
"""
# import neurosynth as ns
# ns.dataset.download(path='.', unpack=True)

from neurosynth.base.dataset import Dataset
from neurosynth import decode

# dataset = Dataset('database.txt')
# dataset.add_features('features.txt')
# dataset.save('dataset.pkl')
dataset = Dataset.load('dataset.pkl')

terms = []
with open("g1_z_pos.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg1_t2_z_cluster_pos.nii', save='real_r_g1_pos.txt')

terms = []
with open("g1_z_neg.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg1_t2_z_cluster_neg.nii', save='real_r_g1_neg.txt')


terms = []
with open("g2_z_pos.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg2_t2_z_cluster_pos.nii', save='real_r_g2_pos.txt')

terms = []
with open("g2_z_neg.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg2_t2_z_cluster_neg.nii', save='real_r_g2_neg.txt')


terms = []
with open("g3_z_pos.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg3_t2_z_cluster_pos.nii', save='real_r_g3_pos.txt')

terms = []
with open("g3_z_neg.txt","r") as f:
    for line in f:
        terms.append(line.strip('\n'))
   
decoder = decode.Decoder(dataset,features=terms)
data = decoder.decode('rg3_t2_z_cluster_neg.nii', save='real_r_g3_neg.txt')
