#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:59:22 2021

@author: feldaher

Aim: Generate random trajectories based experimental data from single neuron tracking
as a "negative control" for trajectory similitude measurements

"""

import os
import itertools
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines


from tkinter import Tk
from tkinter.filedialog import askopenfilename
#from sympy import Point, Line
import numpy as np
import math


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]

#Load experimental data from neurons tracking


# Read track.csv files
dataTrack=pd.read_csv(filenameTracks,usecols=['CID','TID', 'Points'])

# Extract TID and nb points
trackCID=dataTrack["CID"]
tid=dataTrack["TID"]
nPts=dataTrack["Points"]


# Read points.csv files
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenamePts= askopenfilename() # show an "Open" dialog box and return the path to the selected file
truncate_filename=os.path.splitext(os.path.basename(filenamePts))

dataPoints=pd.read_csv(filenamePts,usecols=['CID','x [micron]','y [micron]'])
x=dataPoints["x [micron]"]
y=dataPoints["y [micron]"]
CID=dataPoints["CID"]
#TID=dataPoints["TID"]


nbTraj=dataTrack.shape[0]

# Separate trajectories in the right format for similarity analysis
trajecListNeurons=[]
#distTrajec=pd.DataFrame(index=range(nbTraj), columns=range(nPts))
#distTrajec=pd.DataFrame(index=range(nbTraj), columns=range(nPts[i]))
trajecListSimul=[]

index=0
neuronindex=0


#Initialize random generators
rng = np.random.default_rng()

nbPtsNeuronsSimul=0

for i in range(nbTraj):
	trajec1=pd.DataFrame(data=[dataPoints.loc[index:index+dataTrack.loc[i,"Points"]-1,'x [micron]'],dataPoints.loc[index:index+dataTrack.loc[i,"Points"]-1,'y [micron]']])
	trajec1=trajec1.transpose()
	

	#create list of neurons trajectory
	if dataTrack.loc[i,"CID"]==1:
		trajecListNeurons.append(trajec1.to_numpy())
		
		#Compute distance between successive points of each trajectory
		distance=np.ones(trajec1.shape[0]-1)

		for j in range(trajec1.shape[0]-1):
			nextj=j+1
			distance[j]=math.sqrt((trajec1.iloc[nextj,0]-trajec1.iloc[j,0])**2+(trajec1.iloc[nextj,1]-trajec1.iloc[j,1])**2)
			#print(distance)
#			distTrajec.loc[i,j]=distance

	#Estimate distribution parameters as modeled by a normal distribution
	#Compute average, sigma for distance between points for each track
		print(neuronindex)

		meanDistance=np.mean(distance)
		stdDistance=np.std(distance)
		
		trajecSimul=trajec1.copy()
			
	#Generate array of random values following average and sigma from experimental data. 
	# USE experimental starting points

		for k in range(trajecSimul.shape[0]-1):
			nbPtsNeuronsSimul+=1

			newDistance = np.random.default_rng().normal(meanDistance, stdDistance)
			randAngle = rng.integers(low=0, high=360)
			newX=newDistance*math.cos(math.radians(randAngle))
			newY=newDistance*math.sin(math.radians(randAngle))

			trajecSimul.iloc[k+1,0]=trajecSimul.iloc[k,0]+newX #start at k=1 to keep the same starting points than the original trajectory
			trajecSimul.iloc[k+1,1]=trajecSimul.iloc[k,1]+newY
		
		
		trajecListSimul.append(trajecSimul.to_numpy())
#		neuronCID[neuronindex]=1
#		fig, ax = plt.subplots()
#		q = ax.plot([row[0] for row in trajecListNeurons[neuronindex]],[row[1] for row in trajecListNeurons[neuronindex]])
#		r=ax.plot([row[0] for row in trajecListSimul[neuronindex]],[row[1] for row in trajecListSimul[neuronindex]])

#		plt.show()

		neuronindex+=1

	index=index+dataTrack.loc[i,"Points"]

count = sum( [ len(listElem) for listElem in trajecListSimul])
print(count)

nbPtsNeuronsSimul=nbPtsNeuronsSimul+neuronindex

#Save simulated trajectories for further analysis
#Use file format from MTrackJ

df = pd.DataFrame(np.concatenate(trajecListSimul),columns=['x [micron]','y [micron]'])
df['CID']=np.ones(nbPtsNeuronsSimul).transpose()

df.to_csv(directory+'/'+truncate_filename[0]+'simulated.csv', index = False, header=True)
