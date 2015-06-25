#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is the "snake" Python module. It is the lower level module used to
analyze non-magnetic bright points.

It is assume that small UIO 3D boxes containing single nMBPs are given.
The goal of this module is to locate the 3D region in which nMBP forms.

The set of minima of density (layer by layer, at constant z) correponding
to the nMBP are first located. This set has the shape of a "snake", and
hence the name of the module. This "snake" defines the "eye" of the
nMBP, looking like the "eye" of a tornado.

At each layer, one has the "eye" of the nMBP. For a specific layer and from
the corresponding eye one then computes the region that contains the nMBP.

First step
----------
First step is to obtain candidates for the local minimum (of density), that
corresponds to the eye of the nMBP in a specific layer. For this purpose,
a minimum filter is used.

This introduces a free parameter, 'radius', or 'r'. Experimentally 3 is a
good value.

Unfortunately this gives many candidates and only one has to be selected.

Second step
-----------
The second step is to reduce the number of selected candidates.

The eye of the nMBP has a quite low density, but only around tau=1. For a
layer around tau=1, one can plot a histogram of density at local minima.
Experimentally one obtains that the steepest slope of the histogram
separates fairly well a small set of local minima (among which the nMBP)
from a set containing the biggest amount of irrelevant local minima.

This procedure works extremely well at tau=1, usually leading directly to
the single nMBP. Away from tau=1 it generally leaves a few other
candidates.

The candidates left are sorted by increasing density, and at most the first
three are kept.

Third step
----------
Third step is identifying THE nMBP among all remaining local minima. This
cannot be done generally from one single layer. The idea is to first
apply step two to a layer near tau=1 and for each candidate move up and
down to the nearest local minimum, thus identifying all possible "snakes".

At this point, one has as many snakes as initial candidates. One chooses
the longest, and if two snakes have identical length, one chooses the one
that is less twisted.

Getting the neighbourhood
-------------------------

One needs now to identify the neighbourhood of the nMBP. We proceed again
layer by layer. We look for a star domain whose centre is the already
computed eye.

The idea is to move radially in each direction away from the eye. Density
generally increases, and one always keep track of the maximum density.
When the density increases above a fraction (1.-threshold) of
rho_max - rho_eye, one decides it is the boundary of the star domain.

The threshold is experimentally taken to be 0.8. It is the other free
parameter of the algorithms used here. Note that the results are not
very sensitive to these free parameters.

Of course, the number of directions to explore radially increases with
the size of the domain, and for points near the eye there is a huge
redundancy.

The idea is to build a (discrete) spiral around the eye, parametrized so
that all points are more or less equidistant from the neighbours and that
the maximum distance does not exceed a fraction of the inter-cell
distance. This fraction is chosen so that the spiral covers all the mesh.

This spiral is an invertible mathematical function, and one can therefore
get for each point the neighbours that are in a radial direction,
allowing hence to move radially.

Because the mesh is discrete, the algorithm does not work perfectly
and the star domain found is in fact not even connected, but it seems
to be a very good choice of neighbourhood.

Getting the region containing the nMBP, its border, and the background
----------------------------------------------------------------------

Now one looks for a domain inside the neighbourhood, containing the
nMBP, and satisfying condition C defined by

              |(rho_eye-<rho>)/<rho>| < 0.

This domain must be as big as possible satisfying these conditions.
The background is defined as the complementary part of this region with
respect to the neighbourhood.

This domain D is choosen consistently so that <rho> is the mean taken
over the background and such that the domain is connected. One proceeds
iterarively. A first approximation to D is given by the neighbourhood
itself.

A new approximation is computed applying D followed by a flood fill
algorithm (flood fill produces the biggest connected domain to the eye
and satisfying some condition).

Applying this condition over and over reduces the size of D. When the
size remains constant the algorithm is stopped.
"""

import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
try:
	import matplotlib.pyplot as plt
except ImportError: print('Warning: matplotlib is missing.')
from flood_fill import flood_fill

def find_min(data, ng_sz=3, threshold=np.Inf, time=-1):
	data_max = filters.maximum_filter(data, ng_sz)
	maxima = (data == data_max)
	data_min = filters.minimum_filter(data, ng_sz)
	minima = (data == data_min)
	diff = ((data_max - data_min)/np.abs(data_max) < threshold)
	minima[diff == False] = False
	labeled, num_objects = ndimage.label(minima)
	slices = ndimage.find_objects(labeled)
	#x, y = [], []
	p = []
	for dy, dx in slices:
		x_center = .5*(dx.start+dx.stop-1)
		#x.append(x_center)
		y_center = .5*(dy.start+dy.stop-1)
		#y.append(y_center)
		if time >= 0: p.append((y_center,x_center,time))
		else: p.append((y_center,x_center))

	return p

def x_coords(p):
	return [p0[0] for p0 in p]

def y_coords(p):
	return [p0[1] for p0 in p]

def d2(p1,p2,periodX=0,periodY=0):
	if periodX == 0: deltaX = p2[0]-p1[0]
	else : deltaX = min(np.abs(p2[0]-p1[0]),periodX-np.abs(p2[0]-p1[0])+1)
	if periodY == 0: deltaY = p2[1]-p1[1]
	else : deltaY = min(np.abs(p2[1]-p1[1]),periodY-np.abs(p2[1]-p1[1])+1)
	return (deltaX)**2+(deltaY)**2

def copy_snake(p):
	return [tuple(p0) for p0 in p]

def copy_all_snakes(s):
	return [copy_snake(s0) for s0 in s]

def split_data(v):
	#if len(v) == 0: return 0
	#if len(v) == 1: return v[0]
	#v0_arg = np.argsort(v)
	#v0 = np.array(v)[v0_arg]
	#d = v0[1:]-v0[:-1]
	#idx = np.argsort(d)[-1]
	#return .5*(v0[idx+1]+v0[idx])

	if len(v) == 0: return 0
	h, q = np.histogram(v)
	d = h[1:]-h[:-1]
	idx = np.argsort(d)[-1]
	return .5*(q[idx+1]+q[idx])

def order_snake(p, data):
	#median = filters.median_filter(data,size=10)
	#p_img = [-data[p0[0],p0[1]]/median[p0[0],p0[1]] for p0 in p]
	#mean = filters.generic_filter(data,np.mean,size=10)
	#p_img = [-data[p0[0],p0[1]]/mean[p0[0],p0[1]] for p0 in p]
	p_img = np.array([data[p0[0],p0[1]] for p0 in p])
	keep = np.where(p_img<split_data(p_img))
	p_img = p_img[keep]
	#p_img = [data[p0[0],p0[1]]/np.mean(data[select_disk(data,(int(p0[0]),int(p0[1])))[0]]) for p0 in p]
	p_sort = np.array(p)[keep][np.argsort(p_img)]
	if len(p_sort) > 0:
		if len(p_sort[0]) >= 3:
			return [(p0[0], p0[1], int(p0[2])) for p0 in p_sort]
	return [tuple(p0) for p0 in p_sort]

def snake_add_point(s, p, time, periodX=0, periodY=0, radius=1, direction=1, add_new=True):
	s_copy = copy_all_snakes(s)
	s_changed = np.zeros(len(s), dtype=int)
	is_p_new = np.ones(len(p), dtype=int)
	if direction >= 0:
		pidx = -1
		ptime = -1
	else:
		pidx = 0
		ptime = 1
	for j in range(len(p)):
		p0 = p[j]
		for i in range(len(s_copy)):
			s0 = s[i]
			s0c = s_copy[i]
			if d2(s0c[pidx],p0,periodX,periodY) <= radius**2 and s0c[pidx][2]==time+ptime:
				if not s_changed[i]:
					if direction > 0:
						s0.append(p0)
					if direction < 0:
						s0.insert(0,p0)
					s_changed[i] = 1
				#else:
				elif add_new:
					p_new = copy_snake(s0c)
					if direction > 0:
						p_new.append(p0)
						s.append(p_new)
					if direction < 0:
						p_new.insert(0,p0)
						s.insert(0,p_new)
				elif d2(s0c[pidx],p0,periodX,periodY) < d2(s0c[pidx],s0[pidx],periodX,periodY):
					s0.pop(pidx)
					if direction > 0:
						s0.append(p0)
					if direction < 0:
						s0.insert(0,p0)
				is_p_new[j]=0
	new_loc = np.where(is_p_new)
	if len(new_loc) > 0: new_loc = new_loc[0]
	if add_new:
		for j in new_loc: s.append([p[j]])

def snake_twist(p,periodX=0,periodY=0):
	all_twist2 = [d2(p[i],p[i+1],periodX,periodY) for i in range(len(p)-1)]

	return np.sqrt(np.sum(all_twist2))/len(p)

def optimal_snake(s,periodX=0,periodY=0):
	if len(s) == 0: return s
	sz = np.array([len(s0) for s0 in s])
	s_arr = np.array(s)
	big_s = s_arr[np.where(sz==np.max(sz))]
	twist = np.array([snake_twist(p,periodX,periodY) for p in big_s])

	return big_s[np.where(twist==np.min(twist))]

def snake_from_box(box, radius=1, periodic=True, start=-1):
	if periodic :
		pX = np.shape(box)[0]
		pY = np.shape(box)[1]
	else :
		pX = 0
		pY = 0
	snake = []
	if start<0:
		for i in range(np.shape(box)[2]):
			data = box[:,:,i]
			p = find_min(data, time=i)
			snake_add_point(snake, p, i, pX, pY, radius)
	else:
		data = box[:,:,start]
		p = find_min(data, time=start)
		p = order_snake(p,data)
		p = p[:min(3,len(p))]
		snake = [[p0] for p0 in p]
		for i in range(start+1,np.shape(box)[2]):
			data = box[:,:,i]
			p = find_min(data, time=i)
			snake_add_point(snake, p, i, pX, pY, radius, add_new=False)
		for i in range(start-1,-1,-1):
			data = box[:,:,i]
			p = find_min(data, time=i)
			snake_add_point(snake, p, i, pX, pY, radius, direction=-1, add_new=False)
	#return snake
	opt = optimal_snake(snake,pX,pY)
	if len(opt) > 0:
		if len(opt) > 1:
			print("Warning: multiple snakes found!")
		return np.array(opt[0])
	else: return []

def show_min(data, ng_sz=3, threshold=np.Inf, time=-1):
	p = find_min(data, ng_sz=ng_sz, threshold=threshold, time=time)
	values = np.array([data[p0[0],p0[1]] for p0 in p])
	pmin = list(np.array(p)[np.where(values<split_data(values))])
	#print([data[p0[0],p0[1]] for p0 in p])
	#print([data[p0[0],p0[1]] for p0 in pmin])

	#plt.hist([data[p0[0],p0[1]] for p0 in p])
	#plt.show()
	plt.imshow(data, origin='bottom', interpolation='bilinear')
	plt.plot(y_coords(p), x_coords(p), 'ro')
	plt.plot(y_coords(pmin), x_coords(pmin), 'go')

	plt.show()

def test_snake():
	n=100
	s=100
	ng_sz=3
	step=.05
	N=10
	m=5
	alpha_max=2.*np.pi/12.
	threshold=np.Inf

	y_coord = np.repeat(range(n),n).reshape((n,n))
	x_coord = y_coord.T

	gauss2d = lambda x, y : np.exp(-.5*((x_coord-n*x)**2+(y_coord-n*y)**2)/s)
	data3d = np.zeros((n,n,N))
	x0 = np.empty(m)
	y0 = np.empty(m)
	alpha = np.empty(m)
	for i in range(m):
		x0[i] = np.random.rand()
		y0[i] = np.random.rand()
		alpha[i] = 2.*np.pi*np.random.rand()
	for j in range(N):
		for i in range(m):
			alpha[i] += alpha_max*np.random.rand()
			x0[i] = (x0[i]+step*np.cos(alpha[i]))%1
			y0[i] = (y0[i]+step*np.sin(alpha[i]))%1
			data3d[:,:,j] -= gauss2d(x0[i],y0[i])
	#data -= gauss2d(.1,.2)
	#data -= gauss2d(.2,.7)
	#data -= gauss2d(.5,.5)
	#data -= gauss2d(.6,.8)
	#data -= gauss2d(.9,.1)

	snake = []
	for i in range(np.shape(data3d)[2]):
		data = data3d[:,:,i]
		plt.imshow(data, origin='bottom', interpolation='bilinear')
		p = find_min(data, ng_sz, threshold, time=i)
		snake_add_point(snake, p, i, n, n, radius=10)
		plt.plot(y_coords(p), x_coords(p), 'ro')

		plt.show()

def spiral(t, c0, s):
	c = .3
	r = np.sqrt(c*t/np.pi)
	a = 2.*np.sqrt(np.pi*c*t)
	x = int(round(r*np.cos(a) + c0[1]))
	y = int(round(r*np.sin(a) + c0[0]))
	#x = r*np.cos(a) + c0[1]
	#y = r*np.sin(a) + c0[0]
	if y >= s[0] or x >= s[1]: return None
	return (y,x)

def prev_point(p, c0, s):
	c = .3
	(y,x) = p
	r = np.sqrt((x-c0[1])**2+(y-c0[0])**2)-1.
	deltaX = float(x)-c0[1]
	deltaY = float(y)-c0[0]
	#r = np.sqrt(d2(p,c0,periodX=s[0],periodY=s[1]))-1.
	#deltaX = (float(x)-c0[1]+s[1]/2.)%s[1]-s[1]/2.
	#deltaY = (float(y)-c0[0]+s[0]/2.)%s[0]-s[0]/2.
	if x-c0[1] != 0:
		p = deltaY/deltaX
		a0 = np.mod((1-np.sign(deltaX))/2*np.pi+np.arctan(p),2.*np.pi)
	else: a0 = np.sign(deltaY)*np.pi/2.
	a = max(int(r-1./c)*2.*np.pi+a0,0)
	t = int(round((a/2.)**2/(np.pi*c)))
	return spiral(t,c0,s)
	#r = np.sqrt((x-c0[1])**2+(y-c0[0])**2)
	#if x-c0[1] != 0:
	#	p = (float(y)-c0[0])/(float(x)-c0[1])
	#	a = np.mod((1-np.sign(x-c0[1]))/2*np.pi+np.arctan(p),2.*np.pi)
	#else: a = np.sign(y-c0[0])*np.pi/2.
	#new_r = r-1./max(abs(np.sin(a)), abs(np.cos(a)))
	#if new_r < 0. : new_r = 0.
	#return (int(round(c0[0]+new_r*np.sin(a))), int(round(c0[1]+new_r*np.cos(a))))

def star_domain(data, c0, threshold=.8, extrema='min'):
	if extrema == 'max': data0 = -data
	else: data0 = data
	s = np.shape(data0)
	local_min = data0[c0]
	mask = np.zeros_like(data0, dtype=bool)
	mask[c0] = True
	data_max = np.empty_like(data0)
	data_max[c0] = data0[c0]
	N = int(np.pi*(np.min(np.shape(data0))/4.)**2)
	for i in range(N):
		p1=spiral(i,c0,s)
		if not p1: continue
		p0=prev_point(p1,c0,s)
		data_max[p1] = max(data_max[p0],data0[p1])
		if mask[p0] and data0[p1] > threshold*data_max[p1]+(1.-threshold)*local_min:
			mask[p1] = True
	return mask

def test_spiral(data,c0):
	plt.imshow(data,interpolation='none')
	s = np.shape(data)
	#N = int(np.pi*(np.max(np.shape(data))/2.)**2)
	points = []
	N = 300
	for i in range(N):
		p = spiral(i,c0,s)
		points.append(p)
		plt.plot(p[1],p[0],'o')
	points = np.array(points).T
	plt.plot(points[1],points[0])
	plt.show()

def test_spiral2():
	d=np.zeros((100,100),dtype=bool)
	d[50,50]=True
	for i in range(10000):
		p=spiral(i,(50,50),(100,100))
		pp=prev_point(p,(50,50),(100,100))
		d[p]=True
		if not d[pp]: print('Warning: '+str(pp))

def test_star_domain(data,c0):
	mask = star_domain(data,c0)
	plt.imshow(data,interpolation='none')
	plt.imshow(mask,interpolation='none',alpha=0.7)
	plt.plot(c0[1],c0[0],'o')
	plt.show()

def select_disk(data, c0, mask_flag=True, threshold=.5):
	npts = 0
	mesh = np.array(np.where(np.ones_like(data,dtype=bool)))
	pts = [mesh[0], mesh[1]]
	data_min = data[c0[0], c0[1]]
	if mask_flag: mask = star_domain(data,c0)
	else: mask = np.ones_like(data,dtype=bool)
	background = np.where(mask)
	while len(pts[0]) != npts:
		npts = len(pts[0])
		data_mean = np.mean(data[background])
		data_std = np.std(data[background])
		selection = [i for i in range(len(mesh[0])) if data[mesh[0][i],mesh[1][i]]<=threshold*data_mean+(1.-threshold)*data_min and mask[mesh[0][i],mesh[1][i]]]
		selection = tuple(mesh[:,selection])
		selection_mask = np.zeros_like(data,dtype=bool)
		selection_mask[selection] = True
		selection, border = flood_fill(selection_mask, c0, mask)
		background = [i for i in range(len(mesh[0])) if not selection[mesh[0][i],mesh[1][i]] and mask[mesh[0][i],mesh[1][i]]]
		background = tuple(mesh[:,background])
		selection = np.where(selection)

	return selection, border, background
