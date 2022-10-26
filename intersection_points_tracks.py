#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:15:42 2020

@author: feldaher

Aim: compute all the intersection points from processed track data (cf tracking_vector_analysis_3D.py)
Input: file with x,y,u,v coordinates
"""

import numpy as np
import os
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from sympy import Point, Line

import itertools

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]
   

# Read processed.csv files

dataTrack=pd.read_csv(filenameTracks)

x1=dataTrack["xi"].values
y1=dataTrack["yi"].values
u=dataTrack["u"].values
v=dataTrack["v"].values
x2=x1+u
y2=y1+v



intersectionPoints = pd.DataFrame()
						 
lines=[]							 

for i in range(x1.size):
	lines.append(Line(Point(x1[i], y1[i]),Point(x2[i], y2[i])))
	
intersections=[]


for L1,L2 in itertools.permutations(lines,2):
	intersect=L1.intersection(L2)
	if intersect:
		intersections.append(intersect)
		

for i in range(len(intersections)):
	x_val=(intersections[i][0].evalf()).x
	y_val=(intersections[i][0].evalf()).y
	if ((x_val<210)&(x_val>0)&(y_val<210)&(y_val>0)): # removing intersection points outsite of the tissue, set to 210 microns
		inters=pd.DataFrame({'x':[x_val], 'y':[y_val]})
		intersectionPoints=intersectionPoints.append(inters,ignore_index=True)




#Create new Dataframe for storing coordinates of intersection points and save them in a csv file

truncate_filename=os.path.splitext(os.path.basename(filenameTracks))
intersectionPoints.to_csv (directory+'/'+truncate_filename[0]+'_intersections.csv', index = False, header=True)



