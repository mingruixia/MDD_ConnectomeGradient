# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 14:57:34 2022

@author: mingr
"""
# Libraries
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
 
# Set data
df = pd.DataFrame({
'group': ['A'],
'VIS': [0.1247],
'SMN': [0.145],
'DAN': [0.0958],
'VAN': [0.0823],
'LIM': [0.0825],
'FPN': [0.1714],
'DMN': [0.2653],
'SUB': [0.0375]
})
 
# number of variable
categories=list(df)[1:]
N = len(categories)
 
# We are going to plot the first line of the data frame.
# But we need to repeat the first value to close the circular graph:
values=df.loc[0].drop('group').values.flatten().tolist()
values += values[:1]
values
 
# What will be the angle of each axis in the plot? (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]
 
# Initialise the spider plot
plt.figure(dpi=240)
ax = plt.subplot(111, polar=True)
 
# Draw one axe per variable + add labels

plt.xticks(angles[:-1], categories, color='black', size=18)
ax.tick_params(axis='both',which='major',pad=10)
 
# Draw ylabels
ax.set_rlabel_position(90)
plt.yticks([0.1,0.2,0.3], ["10%","20%","30%"], color="grey", size=10)
plt.ylim(0,0.3)
 
# Plot data
ax.plot(angles, values, linewidth=1, linestyle='solid')
 
# Fill area
ax.fill(angles, values, 'b', alpha=0.1)

# Show the graph
plt.show()

