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
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import shapely.geometry as geom
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from sympy import Point, Line
import numpy as np
import mpl_scatter_density
import math
 

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]

threshold_vector=2.0 # keep only neuron with a significant displacement

ImageDim_x=194.15
ImageDim_y=186.5
intertectal_angle=0.21
x_ref=28.2
y_ref=172.8
x_c=81.7
y_c=127.8


   
'''
@@@@@@@@@@@@@@@@@   Preprocess MTrackJ exported files
'''

# Read track.csv files

dataTrack=pd.read_csv(filenameTracks,usecols=['TID', 'Points'])



# Extract TID and nb points
tid=dataTrack["TID"]
nPts=dataTrack["Points"]


# Read points.csv files
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenamePts= askopenfilename() # show an "Open" dialog box and return the path to the selected file



dataPoints=pd.read_csv(filenamePts,usecols=['Nr','CID','x [µm]','y [µm]','z [µm]'])
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
	if i>0:# Filtering the track with enough points
		if dataPoints["CID"][index]>0: # Filtering tracks to only keep the ones near the injury (CID=2)
			#read and store xyz init and final values
			coordinates.loc[index,"xi"]=x.loc[index+int((i-1)/3)] #Take mid-point as starting point of vectors
			coordinates.loc[index,"yi"]=y.loc[index+int((i-1)/3)]
			coordinates.loc[index,"zi"]=z.loc[index+int((i-1)/3)]
			coordinates.loc[index,"xf"]=x.loc[index+i-1]
			coordinates.loc[index,"yf"]=y.loc[index+i-1]
			coordinates.loc[index,"zf"]=z.loc[index+i-1]
			# Compute vectors
			trackDirection.loc[index,"u"]=coordinates.loc[index,"xf"]-coordinates.loc[index,"xi"]
			trackDirection.loc[index,"v"]=coordinates.loc[index,"yf"]-coordinates.loc[index,"yi"]
			trackDirection.loc[index,"w"]=coordinates.loc[index,"zf"]-coordinates.loc[index,"zi"]
	
		index=index+i


combinedData=pd.concat([coordinates["xi"],coordinates["yi"],coordinates["zi"], trackDirection], axis=1)	

'''
@@@@@@@@@@@@@@@@@   Compute  vectors for displacement field analysis
'''
# Filter vector amplitude to keep only cells with a significant displacement in x or y

filteredData=combinedData[combinedData.u.ge(threshold_vector)|combinedData.v.ge(threshold_vector)|combinedData.u.le(-threshold_vector)|combinedData.v.le(-threshold_vector)]


x1=pd.to_numeric(filteredData["xi"])
y1=pd.to_numeric(filteredData["yi"])
u=pd.to_numeric(filteredData["u"]) 
v=pd.to_numeric(filteredData["v"]) 




Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameSkin= askopenfilename() # show an "Open" dialog box and return the path to the selected file

   
#read skin outline coordinates from ImageJ manual freehand ROI
skin_coords = np.loadtxt(filenameSkin)



# Flip vertical
y1=ImageDim_y-y1
y_ref=ImageDim_y-y_ref
skin_coords[:,1]=ImageDim_y-skin_coords[:,1]
#intertectal_angle=-intertectal_angle
v=-v

#Rotate to orientate the frame of reference verticaly on LSM880 images

x0=ImageDim_x/2
y0=ImageDim_y/2


x2=(y1-y0)+x0
y2=-(x1-x0)+y0

x_ref2=(y_ref-y0)+x0
y_ref2=-(x_ref-x0)+y0


#Rotate vectors
u2=v
v2=-u



#Rotate skin outline
new_skin_coords=skin_coords.copy()
new_skin_coords[:,0]=(skin_coords[:,1]-y0)+x0
new_skin_coords[:,1]=-(skin_coords[:,0]-x0)+y0

#intertectal_angle=intertectal_angle-90




#Flip vertical
y2=ImageDim_x-y2
y_ref2=ImageDim_x-y_ref2
new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]
#intertectal_angle=-intertectal_angle
v2=-v2

# y_ref2=ImageDim_x-y_ref2
# new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]


# Convert into radians
intertectal_angle=math.radians(-intertectal_angle)




fig, ax = plt.subplots()
q = ax.quiver(x2,y2,u2,v2,units='width')
plt.show()


# Write results in file
truncate_filename=os.path.splitext(os.path.basename(filenamePts))
filteredData.to_csv (directory+'/'+truncate_filename[0]+'_processed.csv', index = False, header=True)
fig.savefig(directory+'/'+truncate_filename[0]+'_quiver.svg', format='svg', dpi=600)


'''
@@@@@@@@@@@@@@@@@   Compute  intersection points 
'''
x3=x2+u2
y3=y2+v2

x2=x2.values
y2=y2.values
x3=x3.values
y3=y3.values

intersectionPoints = pd.DataFrame()
						 
lines=[]							 

for i in range(x1.size):
	lines.append(Line(Point(x2[i], y2[i]),Point(x3[i], y3[i])))
	
intersections=[]


for L1,L2 in itertools.permutations(lines,2):
	intersect=L1.intersection(L2)
	if intersect:
		intersections.append(intersect)
		

for i in range(len(intersections)):
	x_val=(intersections[i][0].evalf()).x
	y_val=(intersections[i][0].evalf()).y
#	if ((x_val<210)&(x_val>0)&(y_val<210)&(y_val>0)): # removing intersection points outsite of the tissue, set to 210 microns
	inters=pd.DataFrame({'x':[x_val], 'y':[y_val]})
	intersectionPoints=intersectionPoints.append(inters,ignore_index=True)

intersectionPoints.to_csv (directory+'/'+truncate_filename[0]+'_intersections.csv', index = False, header=True)


x_inter=intersectionPoints['x']
y_inter=intersectionPoints['y']
cgx = np.sum(x_inter)/len(x_inter)
cgy = np.sum(y_inter)/len(y_inter)

'''
@@@@@@@@@@@@@@@@@   Compute  CoG and graph intersection points + CoG
'''
#Save the coordinates of the center of gravity

d = {'xCoG': [cgx], 'yCoG': [cgy]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG.csv', index = False, header=True)

#Display scatter plot and save it as SVG file
plt.scatter(x_inter,y_inter);
plt.scatter(cgx, cgy, color='k', marker='+', s=1e4);

plt.savefig(directory+'/'+truncate_filename[0]+'_intersection_CoG.svg', format='svg', dpi=600)

plt.show()

'''
@@@@@@@@@@@@@@@@@   Compute distance and angle CoG - skin intercept 
'''


line = geom.LineString(skin_coords)


point_CoG = geom.Point(cgx, cgy)

# Note that "line.distance(point_CoG)" would be identical
distance=point_CoG.distance(line)

point_on_line = line.interpolate(line.project(point_CoG))
x_skin=point_on_line.x
y_skin=point_on_line.y

cgx=float(cgx)
cgy=float(cgy)

# !! angle inverted on purpose angle = - arctan
angle=np.degrees(np.arctan2((cgy-point_on_line.y),(cgx-point_on_line.x)))

d = {'x_skin': [x_skin], 'y_skin': [y_skin],'distance_CoG_skin': [distance],'angle':[angle]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG_analysis.csv', index = False, header=True)


'''
@@@@@@@@@@@@@@@@@   Transform coordinates of CoG to be at the same origin (intertecal junction caudal point), with the skin perpendicular to the x axis 
'''

# Translation and rotation 
cgx2=(cgx-x_ref)*math.cos(intertectal_angle)-(cgy-y_ref)*math.sin(intertectal_angle)
cgy2=(cgx-x_ref)*math.sin(intertectal_angle)+(cgy-y_ref)*math.cos(intertectal_angle)



# Convert to polar coordinates

r_cog=math.sqrt(cgx2**2+cgy2**2)
theta_cog=np.degrees(np.arctan2(cgy2,cgx2))



# Save the coordinates of the tranformed center of gravity
d = {'xCoG2': [cgx2], 'yCoG2': [cgy2]}
df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG2.csv', index = False, header=True)

d = {'r_cog': [r_cog], 'theta_cog': [theta_cog]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'__CoG2_polar.csv', index = False, header=True)


'''
@@@@@@@@@@@@@@@@@   Transform coordinates of intersection points 
'''

x_inter=intersectionPoints['x']
y_inter=intersectionPoints['y']


x_inter2=(x_inter-x_ref)*math.cos(intertectal_angle)-(y_inter-y_ref)*math.sin(intertectal_angle)
y_inter2=(x_inter-x_ref)*math.sin(intertectal_angle)+(y_inter-y_ref)*math.cos(intertectal_angle)



# Convert to polar coordinates
r=x_inter2.pow(2)+y_inter2.pow(2)
r=r.pow(1./2)

x_inter2=x_inter2.astype(float)
y_inter2=y_inter2.astype(float)

theta =np.degrees(np.arctan2(y_inter2,x_inter2))

# save data to file
d = {'x_inter2': x_inter2, 'y_inter2': y_inter2}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_Intersection_transf.csv', index = False, header=True)

d = {'r': r, 'theta': theta}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_Intersection_transf_polar.csv', index = False, header=True)


'''
@@@@@@@@@@@@@@@@@   Scaling to take into account interindividual variation in tectum size 
'''
# measure distance new origin to skin intersection pt

point_Ref = geom.Point(x_ref, y_ref)

distanceRef=point_Ref.distance(line)


# scale using the distance origin pt-skin intersection pt

x_inter2=x_inter2/distanceRef
y_inter2=y_inter2/distanceRef



d = {'x_inter2': x_inter2, 'y_inter2': y_inter2}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_Intersection_transf_Scaled.csv', index = False, header=True)


r=x_inter2.pow(2)+y_inter2.pow(2)
r=r.pow(1./2)
d = {'r': r, 'theta': theta}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_Intersection_transf_polar_scaled.csv', index = False, header=True)


# Scale center of gravity
cgx2=cgx2/distanceRef
cgy2=cgy2/distanceRef



# Save the coordinates of the tranformed center of gravity
d = {'xCoG2': [cgx2], 'yCoG2': [cgy2]}
df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG2_scaled.csv', index = False, header=True)

r_cog=math.sqrt(cgx2**2+cgy2**2)

d = {'r_cog': [r_cog], 'theta_cog': [theta_cog]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'__CoG2_polar_scaled.csv', index = False, header=True)



