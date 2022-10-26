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
 

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameTracks = askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenameTracks)[0]

threshold_vector=0 # keep only neuron with a significant displacement

ImageDim_x=183.2
ImageDim_y=183.2
intertectal_angle=-5.7
x_ref=39.8
y_ref=135.7
x_c=85
y_c=108.3



   
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
@@@@@@@@@@@@@@@@@   Compute  vectors and starting points
'''
# Filter vector amplitude to keep only cells with a significant displacement in x or y

filteredData=combinedData[combinedData.u.ge(threshold_vector)|combinedData.v.ge(threshold_vector)|combinedData.u.le(-threshold_vector)|combinedData.v.le(-threshold_vector)]


x1=pd.to_numeric(filteredData["xi"])
y1=pd.to_numeric(filteredData["yi"])
u=pd.to_numeric(filteredData["u"]) 
v=pd.to_numeric(filteredData["v"]) 


#Inverse frame of reference since IJ is vertically flipped


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameSkin= askopenfilename() # show an "Open" dialog box and return the path to the selected file

   
#read skin outline coordinates from ImageJ manual freehand ROI
skin_coords = np.loadtxt(filenameSkin)

'''
@@@@@@@@@@@@@@@@@  # Load the Rostrocaudal retinotopic curve
'''


Tk().withdraw() 
filenameRCcurve= askopenfilename() 

   
RCcurve_coords = np.loadtxt(filenameRCcurve)

#invert y coordinates



# Flip vertical
y1=ImageDim_y-y1
y_ref=ImageDim_y-y_ref
skin_coords[:,1]=ImageDim_y-skin_coords[:,1]
y_c=ImageDim_y-y_c
#intertectal_angle=-intertectal_angle
v=-v



RCcurve_coords[:,1]=ImageDim_y-RCcurve_coords[:,1]
x_RCcurve=RCcurve_coords[:,0]
y_RCcurve=RCcurve_coords[:,1]


#Rotate to orientate the frame of reference verticaly on LSM880 images

x0=ImageDim_x/2
y0=ImageDim_y/2


x2=(y1-y0)+x0
y2=-(x1-x0)+y0

x_ref2=(y_ref-y0)+x0
y_ref2=-(x_ref-x0)+y0

x_c1=(y_c-y0)+x0
y_c1=-(x_c-x0)+y0
#Rotate vectors
u2=v
v2=-u



#Rotate skin outline
new_skin_coords=skin_coords.copy()
new_skin_coords[:,0]=(skin_coords[:,1]-y0)+x0
new_skin_coords[:,1]=-(skin_coords[:,0]-x0)+y0

x_RCcurve2=(y_RCcurve-y0)+x0
y_RCcurve2=-(x_RCcurve-x0)+y0



#Flip vertical
y2=ImageDim_x-y2
y_ref2=ImageDim_x-y_ref2
new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]
y_RCcurve2=ImageDim_x-y_RCcurve2
y_c1=ImageDim_x-y_c1
#intertectal_angle=-intertectal_angle
v2=-v2

# y_ref2=ImageDim_x-y_ref2
# new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]


# Convert into radians
intertectal_angle=math.radians(-intertectal_angle)


# # Write results in file

truncate_filename=os.path.splitext(os.path.basename(filenamePts))



'''
@@@@@@@@@@@@@@@@@   Compute size of the tectum (distance REF-skin)
'''

line = geom.LineString(skin_coords)

# measure distance new origin to skin intersection pt

point_Ref = geom.Point(x_ref, y_ref)

distanceRef=point_Ref.distance(line)






'''
@@@@@@@@@@@@@@@@@   Transform coordinates all points (translation and rotation)
'''
x2=x2.to_numpy()
y2=y2.to_numpy()

x3=(x2-x_ref)*math.cos(intertectal_angle)-(y2-y_ref)*math.sin(intertectal_angle)
y3=(x2-x_ref)*math.sin(intertectal_angle)+(y2-y_ref)*math.cos(intertectal_angle)

x_c2=(x_c1-x_ref)*math.cos(intertectal_angle)-(y_c1-y_ref)*math.sin(intertectal_angle)
y_c2=(x_c1-x_ref)*math.sin(intertectal_angle)+(y_c1-y_ref)*math.cos(intertectal_angle)

x_RCcurve3=(x_RCcurve2-x_ref)*math.cos(intertectal_angle)-(y_RCcurve2-y_ref)*math.sin(intertectal_angle)
y_RCcurve3=(x_RCcurve2-x_ref)*math.sin(intertectal_angle)+(y_RCcurve2-y_ref)*math.cos(intertectal_angle)


'''
@@@@@@@@@@@@@@@@@   Transform vectors (rotation only)
'''

u3=u2*math.cos(intertectal_angle)-v2*math.sin(intertectal_angle)
v3=u2*math.sin(intertectal_angle)+v2*math.cos(intertectal_angle)



'''
@@@@@@@@@@@@@@@@@   Scaling to take into account interindividual variation in tectum size 
'''
# scale using the distance origin pt-skin intersection pt

u3=u3/distanceRef
v3=v3/distanceRef

x_RCcurve3=x_RCcurve3/distanceRef
y_RCcurve3=y_RCcurve3/distanceRef

x3=x3/distanceRef
y3=y3/distanceRef

x_c2=x_c2/distanceRef
y_c2=y_c2/distanceRef



'''
@@@@@@@@@@@@@@@@@   Compute projections of the trajectory starting points on the rostrocaudal curve
'''

x_curve=np.empty(x1.size)
y_curve=np.empty(x1.size)

RCcurve_coords[:,0]=x_RCcurve3
RCcurve_coords[:,1]=y_RCcurve3

d = {'x_RCcurve': x_RCcurve3, 'y_RCcurve': y_RCcurve3}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_RCcurve_transf_Scaled.csv', index = False, header=True)


line_RC = geom.LineString(RCcurve_coords)
i=0
for i in range(x1.size):
	point_Track = geom.Point(x3[i], y3[i])
	point_on_line = line_RC.interpolate(line_RC.project(point_Track))
	x_curve[i]=point_on_line.x
	y_curve[i]=point_on_line.y

	


'''
@@@@@@@@@@@@@@@@@   Compute curvilinear position of each neuron on the RC curve by linear interpolation  
'''

position_curve=np.empty(x_curve.size)

i=0
for i in range(x_curve.size):
	point_curve = geom.Point(x_curve[i],y_curve[i])
	position_curve[i] = line_RC.project(point_curve)



'''
@@@@@@@@@@@@@@@@@   Compute local angle of the RC curve and rotate the vectors to linearize the axis  
'''
u3=u3.to_numpy()
v3=v3.to_numpy()

u4=np.empty(u3.size)
v4=np.empty(v3.size)

i=0
for i in range(y_curve.size):
	#get the points coordinates before and after the projected point of each neuron	
	first=np.where(RCcurve_coords[:,1]>y_curve[i])[0]
	first_pt=first[first.size-1]
	last=np.where(RCcurve_coords[:,1]<y_curve[i])[0]
	last_pt=last[0]
	
	#get the angle of the local segment
	dx=RCcurve_coords[last_pt,0]-RCcurve_coords[first_pt,0]
	dy=RCcurve_coords[last_pt,1]-RCcurve_coords[first_pt,1]
	segm_angle=-math.atan2(dy,dx)
		
	# #apply the rotation to have the normals aligned on the y axis
	u4[i]=u3[i]*math.cos(segm_angle)-v3[i]*math.sin(segm_angle)
	v4[i]=u3[i]*math.sin(segm_angle)+v3[i]*math.cos(segm_angle)


'''
@@@@@@@@@@@@@@@@@   Compute local perpendicular lines to the RC curve and their intersection pts  
'''

nbpt_RC=RCcurve_coords[:,0].size-1

x_pt1=np.empty(nbpt_RC)
y_pt1=np.empty(nbpt_RC)
x_pt2=np.empty(nbpt_RC)
y_pt2=np.empty(nbpt_RC)


i=0
for i in range(nbpt_RC):
	# Compute normal vector for each segment
	x1=RCcurve_coords[i,0]
	y1=RCcurve_coords[i,1]
	x2=RCcurve_coords[i+1,0]
	y2=RCcurve_coords[i+1,1]
	dx=x2-x1
	dy=y2-y1
	x_pt1[i]=x1+dx/2.0
	y_pt1[i]=y1+dy/2.0	
	x_pt2[i]=x1+dx/2.0+dy #pt on the normal vector
	y_pt2[i]=y1+dy/2.0-dx




intersectionPoints = pd.DataFrame()
						 
lines=[]		
i=0					 

for i in range(nbpt_RC):
	lines.append(Line(Point(x_pt1[i], y_pt1[i]),Point(x_pt2[i], y_pt2[i])))


	
intersections=[]


for L1,L2 in itertools.permutations(lines,2):
	intersect=L1.intersection(L2)
	if intersect:
		intersections.append(intersect)
		

for i in range(len(intersections)):
	x_val=(intersections[i][0].evalf()).x
	y_val=(intersections[i][0].evalf()).y
	if ((x_val<1)&(x_val>0)&(y_val<1)&(y_val>0)): # removing intersection points outsite of the tissue, set to 210 microns
		inters=pd.DataFrame({'x':[x_val], 'y':[y_val]})
		intersectionPoints=intersectionPoints.append(inters,ignore_index=True)

intersectionPoints.to_csv (directory+'/'+truncate_filename[0]+'_RC_curve_intersections.csv', index = False, header=True)

x_inter=intersectionPoints['x']
y_inter=intersectionPoints['y']

# Convert to polar coordinates
r=x_inter.pow(2)+y_inter.pow(2)
r=r.pow(1./2)

x_inter=x_inter.astype(float)
y_inter=y_inter.astype(float)

theta =np.degrees(np.arctan2(y_inter,x_inter))


d = {'r': r, 'theta': theta}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_RC_curve_intersections_polar.csv', index = False, header=True)




'''
@@@@@@@@@@@@@@@@@   Compute  CoG and graph intersection points + CoG
'''

cgx = np.sum(x_inter)/len(x_inter)
cgy = np.sum(y_inter)/len(y_inter)

#Save the coordinates of the center of gravity

d = {'xCoG': [cgx], 'yCoG': [cgy]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_RCcurve_CoG.csv', index = False, header=True)



# Convert to polar coordinates

r_cog=math.sqrt(cgx**2+cgy**2)
theta_cog=np.degrees(np.arctan2(cgy,cgx))



d = {'r_cog': [r_cog], 'theta_cog': [theta_cog]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_RCcurve_CoG_polar.csv', index = False, header=True)

'''
@@@@@@@@@@@@@@@@@   Saving files and generating figures 
'''
pos_InjuryCenter_curve=line_RC.project(geom.Point(x_c2,y_c2))

curve_length=line_RC.length


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






