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
 




x_ref=59
y_ref=204
ImageDim_y=224.5
ImageDim_x=219
intertectal_angle=-3


   
'''
@@@@@@@@@@@@@@@@@   Preprocess MTrackJ exported files add apply transformations
'''

# Read track.csv files


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenamePts= askopenfilename() # show an "Open" dialog box and return the path to the selected file

directory = os.path.split(filenamePts)[0]
   
#read skin outline coordinates from ImageJ manual freehand ROI
microglia_coords = np.loadtxt(filenamePts)
x=microglia_coords[:,0]
y=microglia_coords[:,1]


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filenameSkin= askopenfilename() # show an "Open" dialog box and return the path to the selected file

   
#read skin outline coordinates from ImageJ manual freehand ROI
skin_coords = np.loadtxt(filenameSkin)



# Flip vertical
y=ImageDim_y-y
y_ref=ImageDim_y-y_ref
skin_coords[:,1]=ImageDim_y-skin_coords[:,1]
#intertectal_angle=-intertectal_angle


#Rotate to orientate the frame of reference verticaly on LSM880 images

x0=ImageDim_x/2
y0=ImageDim_y/2


x1=(y-y0)+x0
y1=-(x-x0)+y0

x_ref2=(y_ref-y0)+x0
y_ref2=-(x_ref-x0)+y0





#Rotate skin outline
new_skin_coords=skin_coords.copy()
new_skin_coords[:,0]=(skin_coords[:,1]-y0)+x0
new_skin_coords[:,1]=-(skin_coords[:,0]-x0)+y0

#intertectal_angle=intertectal_angle-90




#Flip vertical
y1=ImageDim_x-y1
y_ref2=ImageDim_x-y_ref2
new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]
#intertectal_angle=-intertectal_angle


# y_ref2=ImageDim_x-y_ref2
# new_skin_coords[:,1]=ImageDim_x-new_skin_coords[:,1]


# Convert into radians
intertectal_angle=math.radians(-intertectal_angle)




truncate_filename=os.path.splitext(os.path.basename(filenamePts))

'''
@@@@@@@@@@@@@@@@@   Compute  CoG and graph intersection points + CoG
'''

cgx = np.sum(x1)/len(x1)
cgy = np.sum(y1)/len(y1)

#Save the coordinates of the center of gravity

d = {'xCoG': [cgx], 'yCoG': [cgy]}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_CoG.csv', index = False, header=True)

#Display scatter plot and save it as SVG file
plt.scatter(x1,y1);
plt.scatter(cgx, cgy, color='k', marker='+', s=1e4);

plt.savefig(directory+'/'+truncate_filename[0]+'_accumulation_pts.svg', format='svg', dpi=600)

plt.show()

'''
@@@@@@@@@@@@@@@@@   Compute distance and angle CoG - skin intercept 
'''





#print(new_skin_coords[:,0])
line = geom.LineString(new_skin_coords)


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
cgx2=(cgx-x_ref2)*math.cos(intertectal_angle)-(cgy-y_ref2)*math.sin(intertectal_angle)
cgy2=(cgx-x_ref2)*math.sin(intertectal_angle)+(cgy-y_ref2)*math.cos(intertectal_angle)



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
@@@@@@@@@@@@@@@@@   Transform coordinates of accumulation points 
'''


x2=(x1-x_ref2)*math.cos(intertectal_angle)-(y1-y_ref2)*math.sin(intertectal_angle)
y2=(x1-x_ref2)*math.sin(intertectal_angle)+(y1-y_ref2)*math.cos(intertectal_angle)



# Convert to polar coordinates
r=x2**2+y2**2
r=np.sqrt(r)

x2=x2.astype(float)
y2=y2.astype(float)

theta =np.degrees(np.arctan2(y2,x2))

# save data to file
d = {'x_inter2': x2, 'y_inter2': y2}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_accumulation_transf.csv', index = False, header=True)

d = {'r': r, 'theta': theta}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_accumulation_transf_polar.csv', index = False, header=True)


#Display scatter plot and save it as SVG file
plt.scatter(x2,y2);
plt.scatter(cgx, cgy, color='k', marker='+', s=1e4);

plt.savefig(directory+'/'+truncate_filename[0]+'_accumulation_pts_rot.svg', format='svg', dpi=600)

plt.show()

'''
@@@@@@@@@@@@@@@@@   Scaling to take into account interindividual variation in tectum size 
'''
# measure distance new origin to skin intersection pt

point_Ref = geom.Point(x_ref2, y_ref2)

distanceRef=point_Ref.distance(line)
print(distanceRef)


# scale using the distance origin pt-skin intersection pt

x2=x2/distanceRef
y2=y2/distanceRef



d = {'x_inter2': x2, 'y_inter2': y2}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_accumulation_transf_Scaled.csv', index = False, header=True)


r=x2**2+y2**2
r=np.sqrt(r)
d = {'r': r, 'theta': theta}

df = pd.DataFrame(data=d)
df.to_csv(directory+'/'+truncate_filename[0]+'_accumulation_transf_polar_scaled.csv', index = False, header=True)


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



