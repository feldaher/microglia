#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 17:00:34 2020

@author: feldaher

Aim: analysis of intersection points computed by intersection_points_tracks.py
Input: intersection csv file
Output: center of gravity coordinates and scatter plot of intersection points with CoG
"""

import numpy as np
import mpl_scatter_density
import matplotlib.pyplot as plt
import os
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenamePoints= askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenamePoints)[0]
   

# Read processed.csv files

dataTrack = pd.read_csv(filenamePoints)
# Generate fake data


x = dataTrack['x']
y = dataTrack['y']
cgx = np.sum(x)/len(x)
cgy = np.sum(y)/len(y)


#Save the coordinates of the center of gravity
truncate_filename=os.path.splitext(os.path.basename(filenamePoints))

d = {'xCoG': [cgx], 'yCoG': [cgy]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG.csv', index = False, header=True)

#Display scatter plot and save it as SVG file
plt.scatter(x,y);
plt.scatter(cgx, cgy, color='k', marker='+', s=1e4);

plt.savefig(directory+'/'+truncate_filename[0]+'_intersection_CoG.svg', format='svg', dpi=600)

plt.show()

