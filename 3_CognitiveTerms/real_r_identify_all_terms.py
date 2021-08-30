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
decoder = decode.Decoder(dataset)
data = decoder.decode('rg1_t2_z_cluster_pos.nii', save='real_r_g1_pos_all.txt')
data = decoder.decode('rg1_t2_z_cluster_neg.nii', save='real_r_g1_neg_all.txt')
data = decoder.decode('rg2_t2_z_cluster_pos.nii', save='real_r_g2_pos_all.txt')
data = decoder.decode('rg2_t2_z_cluster_neg.nii', save='real_r_g2_neg_all.txt')
data = decoder.decode('rg3_t2_z_cluster_pos.nii', save='real_r_g3_pos_all.txt')
data = decoder.decode('rg3_t2_z_cluster_neg.nii', save='real_r_g3_neg_all.txt')
