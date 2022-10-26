#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:15:49 2020

@author: feldaher
"""

# Aim: extract x,y,z init and final for each track from MTrackJ files
# input: track and points csv files
import os
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from tkinter import Tk
from tkinter.filedialog import askopenfilename


 


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]
   

# Read track.csv files

dataTrack=pd.read_csv(filenameTracks,usecols=['TID', 'Points'])



# Extract TID and nb points
tid=dataTrack["TID"]
nPts=dataTrack["Points"]


# Read points.csv files
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenamePts= askopenfilename() # show an "Open" dialog box and return the path to the selected file

dataPoints=pd.read_csv(filenamePts,usecols=['Nr','TID','x [µm]','y [µm]','z [µm]'])
x=dataPoints["x [µm]"]
y=dataPoints["y [µm]"]
z=dataPoints["z [µm]"]


# Extract for each Track initial and final values of xyz

column_names = ["xi", "yi", "zi","xf", "yf", "zf"]
coordinates= pd.DataFrame(columns = column_names)

column_names_vector = ["u","v","w"]
trackDirection=pd.DataFrame(columns=column_names_vector)


index=0
for i in nPts:
	#read and store xyz init and final values
	coordinates.loc[index,"xi"]=x.loc[index]
	coordinates.loc[index,"yi"]=y.loc[index]
	coordinates.loc[index,"zi"]=z.loc[index]
	coordinates.loc[index,"xf"]=x.loc[index+i-1]
	coordinates.loc[index,"yf"]=y.loc[index+i-1]
	coordinates.loc[index,"zf"]=z.loc[index+i-1]
	# Compute vectors
	trackDirection.loc[index,"u"]=coordinates.loc[index,"xf"]-coordinates.loc[index,"xi"]
	trackDirection.loc[index,"v"]=coordinates.loc[index,"yf"]-coordinates.loc[index,"yi"]
	trackDirection.loc[index,"w"]=coordinates.loc[index,"zf"]-coordinates.loc[index,"zi"]

	index=index+i


 
 
combinedData=pd.concat([coordinates["xi"],coordinates["yi"],coordinates["zi"], trackDirection], axis=1)	



# Create 3d quiver plot and save it
# fig = plt.figure()
# ax = fig.gca(projection='3d')

# ax.quiver(combinedData["xi"], combinedData["yi"],combinedData["zi"], combinedData["u"],combinedData["v"],combinedData["w"], length=0.1, normalize=True)

# Create 2D quiver plot and save it

# Filter vector amplitude to keep only cells with a significant displacement in x or y
threshold_vector=5.0
filteredData=combinedData[combinedData.u.ge(threshold_vector)|combinedData.v.ge(threshold_vector)|combinedData.u.le(-threshold_vector)|combinedData.v.le(-threshold_vector)]


x=pd.to_numeric(filteredData["xi"])
y=pd.to_numeric(filteredData["yi"])
u=pd.to_numeric(filteredData["u"]) 
v=pd.to_numeric(filteredData["v"]) 

#Inverse frame of reference since IJ is vertically flipped
ImageDim_y=209.14
y=ImageDim_y-y
v=-v


fig, ax = plt.subplots()
q = ax.quiver(x,y,u,v,units='width')
plt.show()




# Write results in file
truncate_filename=os.path.splitext(os.path.basename(filenamePts))
filteredData.to_csv (directory+'/'+truncate_filename[0]+'_processed.csv', index = False, header=True)
fig.savefig(directory+'/'+truncate_filename[0]+'_quiver.svg', format='svg', dpi=600)
