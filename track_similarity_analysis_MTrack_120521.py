#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:15:49 2020

@author: feldaher
"""

# Aim: extract x,y,z init and final for each track from MTrackJ files
# input: track and points csv files
import os
import itertools
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import shapely.geometry as geom
import shapely.affinity as affinity
from shapely.ops import nearest_points
from sympy import Point, Line

from tkinter import Tk
from tkinter.filedialog import askopenfilename
#from sympy import Point, Line
import numpy as np
import mpl_scatter_density
import math

#Library for trajectory similarity estimation
import traj_dist.distance as tdist

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]


'''
@@@@@@@@@@@@@@@@@@@ Choice of the metric for estimating trajectory similitude


	1. 'sspd'
        Computes the distances using the Symmetrized Segment Path distance.
    2. 'dtw'
        Computes the distances using the Dynamic Path Warping distance.
    3. 'lcss'
        Computes the distances using the Longuest Common SubSequence distance
    4. 'hausdorff'
        Computes the distances using the Hausdorff distance.
    5. 'frechet'
        Computes the distances using the Frechet distance.
    6. 'discret_frechet'
        Computes the distances using the Discrete Frechet distance.
    7. 'sowd_grid'
        Computes the distances using the Symmetrized One Way Distance.
    8. 'erp'
        Computes the distances using the Edit Distance with real Penalty.
    9. 'edr'
        Computes the distances using the Edit Distance on Real sequence.
		
'''

#similMetric="discret_frechet"

similMetric="hausdorff"


   
'''
@@@@@@@@@@@@@@@@@   Preprocess MTrackJ exported files
'''

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

dataPoints=pd.read_csv(filenamePts,usecols=['Nr','CID','TID','PID','x [micron]','y [micron]'])
x=dataPoints["x [micron]"]
y=dataPoints["y [micron]"]
CID=dataPoints["CID"]
#TID=dataPoints["TID"]


nbTraj=dataTrack.shape[0]

# Separate trajectories in the right format for similarity analysis
trajecListNeurons=[]
trajecListMicroglia=[]

index=0
for i in range(nbTraj):
	trajec1=pd.DataFrame(data=[dataPoints.loc[index:index+dataTrack.loc[i,"Points"]-1,'x [micron]'],dataPoints.loc[index:index+dataTrack.loc[i,"Points"]-1,'y [micron]']])
	trajec1=trajec1.transpose()
	#create list of neurons trajectory
	if dataTrack.loc[i,"CID"]==1:
		trajecListNeurons.append(trajec1.to_numpy())
	#create list for microglia
	elif dataTrack.loc[i,"CID"]==3:
		trajecListMicroglia.append(trajec1.to_numpy())
	
	index=index+dataTrack.loc[i,"Points"]


# Extract for each Track initial and final values of x, y and compute trajectory vectors
column_names = ["xi", "yi","xf","yf"]
neuronCoordinates= pd.DataFrame(columns = column_names)
microgliaCoordinates= pd.DataFrame(columns = column_names)

column_names_vector = ["u","v"]
neuronVector=pd.DataFrame(columns=column_names_vector)
microgliaVector=pd.DataFrame(columns=column_names_vector)



index=0
for i in nPts:
	if dataPoints.loc[index,"CID"]==1:
		neuronCoordinates.loc[index,"xi"]=x.loc[index]
		neuronCoordinates.loc[index,"yi"]=y.loc[index]
		neuronCoordinates.loc[index,"xf"]=x.loc[index+i-1]
		neuronCoordinates.loc[index,"yf"]=y.loc[index+i-1]
		neuronVector.loc[index,"u"]=neuronCoordinates.loc[index,"xf"]-neuronCoordinates.loc[index,"xi"]
		neuronVector.loc[index,"v"]=neuronCoordinates.loc[index,"yf"]-neuronCoordinates.loc[index,"yi"]

	#create list for microglia
	if dataPoints.loc[index,"CID"]==3:
		microgliaCoordinates.loc[index,"xi"]=x.loc[index]
		microgliaCoordinates.loc[index,"yi"]=y.loc[index]
		microgliaCoordinates.loc[index,"xf"]=x.loc[index+i-1]
		microgliaCoordinates.loc[index,"yf"]=y.loc[index+i-1]
		microgliaVector.loc[index,"u"]=microgliaCoordinates.loc[index,"xf"]-microgliaCoordinates.loc[index,"xi"]
		microgliaVector.loc[index,"v"]=microgliaCoordinates.loc[index,"yf"]-microgliaCoordinates.loc[index,"yi"]

	index=index+i


combinedNeuronsCoord=pd.concat([neuronCoordinates["xi"],neuronCoordinates["yi"],neuronVector["u"],neuronVector["v"]], axis=1)	
combinedMicrogliaCoord=pd.concat([microgliaCoordinates["xi"],microgliaCoordinates["yi"],microgliaVector["u"],microgliaVector["v"]], axis=1)	


nbMicroglia=len(combinedMicrogliaCoord.axes[0])
nbNeurons=len(combinedNeuronsCoord.axes[0])

distanceNM = pd.DataFrame(index=range(nbNeurons),columns=range(nbMicroglia))

index=0
for i in range(nbMicroglia):
	for index in range(nbNeurons):
		distanceNM.iloc[index,i]=np.sqrt((combinedMicrogliaCoord.iloc[i,0]-combinedNeuronsCoord.iloc[index,0])**2+(combinedMicrogliaCoord.iloc[i,1]-combinedNeuronsCoord.iloc[index,1])**2)

distanceNM.to_csv(directory+'/'+truncate_filename[0]+'_N_M_distance.csv', index = False, header=True)



''' 
@@@@ Normalise trajectories: translate, rotate, scale to align start and end points 
'''


xi_Neuron=pd.to_numeric(combinedNeuronsCoord["xi"])
yi_Neuron=pd.to_numeric(combinedNeuronsCoord["yi"])
u_Neuron=pd.to_numeric(combinedNeuronsCoord["u"]) 
v_Neuron=pd.to_numeric(combinedNeuronsCoord["v"]) 

xi_Microglia=pd.to_numeric(combinedMicrogliaCoord["xi"])
yi_Microglia=pd.to_numeric(combinedMicrogliaCoord["yi"])
u_Microglia=pd.to_numeric(combinedMicrogliaCoord["u"]) 
v_Microglia=pd.to_numeric(combinedMicrogliaCoord["v"]) 

trajecListNeurons_Transl=[]
trajecListMicroglia_Transl=[]

trajecListNeurons2=trajecListNeurons.copy()
trajecListNeurons3=trajecListNeurons.copy()
	
trajecListMicroglia2=trajecListMicroglia.copy()
trajecListMicroglia3=trajecListMicroglia.copy()

#trajecSimil=pd.DataFrame(index=range(nbMicroglia))

index=0
indexMicroglia=0

#loop on microgliaCoordinates

#initialise trajecListNeurons_Transl

for indexMicroglia in range(nbMicroglia):
	trajecSimil=[]
#	trajecListMicroglia2.clear()
	trajecListNeurons_Transl.clear()
	
	trajecListNeurons2=trajecListNeurons.copy()
	trajecListNeurons3=trajecListNeurons.copy()
	trajecListMicroglia2=trajecListMicroglia.copy()
	trajecListMicroglia3=trajecListMicroglia.copy()

	trajecListMicroglia2[indexMicroglia]=trajecListMicroglia[indexMicroglia]-trajecListMicroglia[indexMicroglia][0]

	
	microgliaAngle=math.atan2(v_Microglia.iloc[indexMicroglia],u_Microglia.iloc[indexMicroglia])

	nbPtsMicroglia=int(trajecListMicroglia[indexMicroglia].size/2.0)
	for j in range(nbPtsMicroglia):
		trajecListMicroglia3[indexMicroglia][j][0]=trajecListMicroglia2[indexMicroglia][j][0]*math.cos(-microgliaAngle)-trajecListMicroglia2[indexMicroglia][j][1]*math.sin(-microgliaAngle)
		trajecListMicroglia3[indexMicroglia][j][1]=trajecListMicroglia2[indexMicroglia][j][0]*math.sin(-microgliaAngle)+trajecListMicroglia2[indexMicroglia][j][1]*math.cos(-microgliaAngle)

for index in range(nbNeurons):
		#Translation - set origin to 0,0 (ref image)
	trajecListNeurons2[index]=trajecListNeurons[index]-trajecListNeurons[index][0]
		#	print("index:" + str(index))
		#	print(*trajecListNeurons[index])
		
		#Rotation
	neuronAngle=math.atan2(v_Neuron.iloc[index],u_Neuron.iloc[index])

	nbPointsNeurons=int(trajecListNeurons[index].size/2.0)
	for k in range(nbPointsNeurons):
		trajecListNeurons3[index][k][0]=trajecListNeurons2[index][k][0]*math.cos(-neuronAngle)-trajecListNeurons2[index][k][1]*math.sin(-neuronAngle)
		trajecListNeurons3[index][k][1]=trajecListNeurons2[index][k][0]*math.sin(-neuronAngle)+trajecListNeurons2[index][k][1]*math.cos(-neuronAngle)
			#print("rotation neurons")

#		scaleFactor=abs(trajecListMicroglia3[indexMicroglia][nbPtsMicroglia-1][0]/trajecListNeurons3[index][nbPointsNeurons-1][0])
#			np.sqrt((u_Microglia.iloc[indexMicroglia]**2+v_Microglia.iloc[indexMicroglia]**2)/(u_Neuron.iloc[index]**2+v_Neuron.iloc[index]**2))
#			trajecListNeurons3[index][j][0]=trajecListNeurons3[index][j][0]*scaleFactor
			
#		trajecListNeurons_Transl.append(trajecListNeurons3[index])                                                                                    
	

	#trajecListMicroglia_Transl.append(trajecListMicroglia2[indexMicroglia])

trajecListNeurons4=trajecListNeurons3.copy()

for indexMicroglia in range(nbMicroglia):
	trajecSimil=[]
	trajecListNeurons4=trajecListNeurons3.copy()
	
	trajecListMicroglia_Transl.clear()
	trajecListNeurons_Transl.clear()

	
	nbPtsMicroglia=int(trajecListMicroglia[indexMicroglia].size/2.0)
	
	for index in range(nbNeurons):
#		print(index)
		nbPointsNeurons=int(trajecListNeurons[index].size/2.0)
		scaleFactor=abs(trajecListMicroglia3[indexMicroglia][nbPtsMicroglia-1][0]/trajecListNeurons3[index][nbPointsNeurons-1][0])
		
		
		for k in range(nbPointsNeurons):
			trajecListNeurons4[index][k][0]=trajecListNeurons3[index][k][0]*scaleFactor
			
		trajecListNeurons_Transl.append(trajecListNeurons4[index])   
	
	trajecListMicroglia_Transl.append(trajecListMicroglia3[indexMicroglia])	

	# COMPUTE TRAJECTORY DISTANCE using traj_dist library (https://github.com/bguillouet/traj-dist/blob/master/traj_dist/distance.py)
	trajecSimil=tdist.cdist(trajecListNeurons_Transl,trajecListMicroglia_Transl,metric=similMetric, type_d="euclidean", converted=None,
          precision=None, eps=None, g=None, verbose=True)

	# SAVE DATA

	df = pd.DataFrame(data=trajecSimil)
	df.to_csv(directory+'/'+truncate_filename[0]+'_'+'microglia'+str(indexMicroglia)+'_'+similMetric+'_Trajectory_distance.csv', index = False, header=True)

	del df




	
	fig, ax = plt.subplots()
	q = ax.plot([row[0] for row in trajecListNeurons4[0]],[row[1] for row in trajecListNeurons4[0]],[row[0] for row in trajecListMicroglia3[indexMicroglia]],[row[1] for row in trajecListMicroglia3[indexMicroglia]])

	# l2 = lines.Line2D([pos_InjuryCenter_curve-0.1, pos_InjuryCenter_curve-0.1], [1.0, 1.01],color='g')
	# ax.add_line(l2)
	plt.show()




'''----- EDIT ---------------------------
create a plot for each microglia


#Center RC curve to injury center
position_curve=position_curve-pos_InjuryCenter_curve
fig, ax = plt.subplots()
q = ax.quiver(position_curve,1,u4,v4,units='width')

# l2 = lines.Line2D([pos_InjuryCenter_curve-0.1, pos_InjuryCenter_curve-0.1], [1.0, 1.01],color='g')
# ax.add_line(l2)
plt.show()

fig.savefig(directory+'/'+truncate_filename[0]+'_curv_quiver_centred.svg', format='svg', dpi=600)

d = {'curvilinear pos':position_curve,'u': u4,'v':v4,'Inj. center curv. coord.':pos_InjuryCenter_curve}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_curve_transf_Scaled_centred.csv', index = False, header=True)


'''








