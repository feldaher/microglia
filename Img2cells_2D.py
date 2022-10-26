#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 17:39:59 2019

@author: fwaharte
"""

#open 3D stack from .tiff multipage file

from skimage.io import imread, imsave
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection
from collections import deque


def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83

    Make a scatter plot of circles.
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_)
               for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    #ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    return collection



img=imread('microglia_cloud2.tif')

print(img.shape)
#plt.imshow(img[0,:,:])

xmin = -500.
xmax = 500
ymin = -500
ymax = 500
x_range = xmax - xmin
y_range = ymax - ymin

yval = 0.0
num_cells = 0

fig_height = 7.0
fig_width = fig_height*x_range/y_range
fig = plt.figure(figsize=(fig_width,fig_height))
ax = fig.gca()
#ax.set_aspect("equal")

xlist = deque()
ylist = deque()
rlist = deque()
rgb_list = deque()

cell_radius = 8
cell_spacing = 1.2 * 2.0 * cell_radius;
print("cell radius and spacing = ",cell_radius,cell_spacing)





img_xrange = img.shape[0]
img_yrange = img.shape[1]


fp = open('cells.dat','w')



ydel = cell_spacing * np.sqrt(3.0)/2.0


#ydel = cell_radius * np.sqrt(3.0)/2.0



ny = 0
zval = 0.0
# Loop over entire (mesh) domain and create hex-packed cells only where image appears.
for yval in np.arange(ymin,ymax,ydel):
	idy = int(img_yrange * (yval-ymin)/y_range)
	xdel = 0.0
	if ny % 2:
		xdel = 0.5*cell_spacing
	ny += 1
	for xval in np.arange(xmin+xdel,xmax+xdel,cell_spacing):
		idx = int(img_xrange * (xval-xmin)/x_range)
		if idx > img_xrange:
			continue
		if (img[idx,idy] >=1):  # red-ish
			xlist.append(yval)
			ylist.append(-xval)
			rlist.append(cell_radius)
			rgb = list(map(int, "255,0,0".split(",")))
			rgb[:]=[x/255. for x in rgb]
			rgb_list.append(rgb)

			cell_type = 1
			cell_str = '%f %f\n' % (yval , -xval)
			fp.write(cell_str)
			num_cells += 1

fp.close()


xvals = np.array(xlist)
yvals = np.array(ylist)
rvals = np.array(rlist)
rgbs = np.array(rgb_list)

#  print('type(rgbs) = ',type(rgbs))
#  print('rgbs = ',rgbs)
#print("xvals[0:5]=",xvals[0:5])
#print("rvals[0:5]=",rvals[0:5])
#  print("rvals.min, max=",rvals.min(),rvals.max())

# plt.cla()
#   title_str += " (" + str(num_cells) + " agents)"
#   plt.title(title_str)
axes_min = xmin
axes_max = xmax
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
circles(xvals,yvals, s=rvals, color=rgbs, edgecolor='black')    # alpha=1.0, edgecolor='black'

plt.show()
