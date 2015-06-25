#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from flood_fill import flood_fill
from snake import find_min
try:
	from matplotlib import pyplot as plt
	from matplotlib import cm
except ImportError: print('Warning: matplotlib is missing.')
try:
	import scipy.ndimage.filters as filters
except ImportError: print('Warning: scipy is missing.')

def getThreshold(data):
	h = np.histogram(data, bins=100)
	return h[1][int(np.argmax(h[0]))]

def plotThreshold(data):
	t = getThreshold(data)
	f, (ax1, ax2) = plt.subplots(1,2)
	ax1.imshow(data, cmap=cm.gray)
	ax2.imshow(data, cmap=cm.gray)
	ax2.imshow(data>t, alpha=.6)
	plt.show()

def circleFootprint(N):
	fp = np.zeros((N,N), dtype=bool)
	r = (N-1.)/2.
	for i in range(N):
		for j in range(N):
			if (i-r)**2+(j-r)**2 <= r**2: fp[i,j] = True
	return fp

def getGranules(data, shrink=.05, verbose=False):
	print('Get Threshold')
	t = getThreshold(data)
	print('np.where')
	granules = np.where(data>t)
	ff = float(np.size(data)-len(granules[0]))/float(np.size(data))
	sz = np.sqrt(np.prod(np.shape(data)))
	fp = circleFootprint(int(np.round(shrink*ff*sz)))
	print('np.zeros_like')
	data0 = np.zeros_like(data, dtype=bool)
	print('fill with True')
	data0[granules] = True
	print('min filter')
	data0 = filters.minimum_filter(data0,footprint=fp)
	print('logical not')
	intergranules = np.logical_not(data0)
	print('np.where')
	ig_pts = np.where(intergranules)
	ff = 0.
	print('select intergranule')
	while ff < .1:
		idx = int(len(ig_pts[0])*np.random.rand())
		print('flood_fill')
		ig, b = flood_fill(intergranules, (ig_pts[0][idx], ig_pts[1][idx]))
		ff = float(len(np.where(ig)[0]))/float(np.size(data))
	print('logical not')
	data0 = np.logical_not(ig)
	print('np.where')
	granules = np.where(data0)
	granules_list = []
	while len(granules[0]) > 0:
		granule, b = flood_fill(data0, (granules[0][0], granules[1][0]))
		data0[granule] = False
		granules = np.where(data0)
		granule_pts = np.where(granule)
		cX = int(round(np.mean(granule_pts[0])))
		cY = int(round(np.mean(granule_pts[1])))
		granules_list += [((cY, cX),granule_pts)]
	granules_sz = [len(g[1][0]) for g in granules_list]
	idx = np.where(granules_sz>(shrink*ff*sz)**2)
	#plt.hist(granules_sz, bins=30)
	granules_list = [granules_list[i] for i in idx[0]]
	if verbose:
		return granules_list
	else:
		granules = np.zeros_like(data, dtype=bool)
		for g in granules_list:
			granules[g[1]] = True
		return granules

def plotGranules(data):
	granules = getGranules(data, verbose=True)
	f, (ax1, ax2) = plt.subplots(1,2)
	ax1.imshow(data, cmap=cm.gray, origin='bottom')
	ax2.imshow(data, cmap=cm.gray, origin='bottom')
	cX = []
	cY = []
	granules_img = np.zeros_like(data, dtype=bool)
	for g in granules:
		granules_img[g[1]] = True
		cX += [g[0][0]]
		cY += [g[0][1]]
		#granule = np.zeros_like(data, dtype=bool)
		#granule[g[1]] = True
		#ax2.imshow(granule, alpha=.1)
		#ax2.plot(g[0][0],g[0][1],'o')
	ax2.plot(cX, cY,'o')
	ax2.imshow(granules_img, alpha=.6, origin='bottom')
	ax2.set_xlim(0, np.size(data, axis=1))
	ax2.set_ylim(0, np.size(data, axis=0))
	plt.show()

def checkNMBP(data, p, granules, size=None, footprint=None):
	if size != None: fp = circleFootprint(size)
	else: fp = footprint
	xshift = int(p[0]-(np.size(fp, axis=0)-1)/2.)
	yshift = int(p[1]-(np.size(fp, axis=1)-1)/2.)
	NX, NY = np.shape(data)
	n = np.count_nonzero(fp)
	mask = np.array(np.where(fp))+np.repeat([[xshift],[yshift]], n, axis=1)
	mask = (np.mod(mask[0], NX), np.mod(mask[1], NY))
	nh = data[mask]*(1-granules[mask].astype(int))	# Neighbourhood
	dnh = nh[np.where(nh<np.mean(nh))] # Dark neighbourhood
	#if np.mean(data[mask])/data[p] > .8 and not np.any(granules[mask]):
	#S = float(np.count_nonzero(fp))
	#if np.mean(dnh)/data[p] < .75 and np.count_nonzero(granules[mask])/S < .1:
	if np.mean(dnh)/data[p] < .85 and not np.any(granules[mask]):
		return True
	else: return False
	#return mask

def getNMBPs(data, granules=None):
	if granules == None: granules=getGranules(data)
	p = find_min(-data)
	p = [p0 for p0 in p if checkNMBP(data, p0, granules, 20)]
	#p = [p0 for p0 in p if not granules[p0]]
	print(p)
	return p

def showNMBPs(data, granules=None):
	if granules == None: granules=getGranules(data)
	p = getNMBPs(data, granules)
	ig_int = np.mean(data[np.where(np.logical_not(granules))])
	plt.hist([data[p0] for p0 in p]/ig_int, bins=30)
	h = np.histogram([data[p0] for p0 in p]/ig_int, bins=30)
	threshold = h[1][np.argmin(h[0])]
	#p = [p0 for p0 in p if data[p0]/ig_int>threshold]
	p = [p0 for p0 in p if data[p0]/ig_int>1.]
	fig, (ax1, ax2) = plt.subplots(1, 2)
	ax2.imshow(data,origin='bottom',cmap=cm.gray)
	ax2.imshow(granules,origin='bottom',cmap=cm.gray,alpha=.2)
	ax2.plot([p0[1] for p0 in p],[p0[0] for p0 in p], 'o')
	ax2.set_xlim(0, np.size(data, axis=1))
	ax2.set_ylim(0, np.size(data, axis=0))
	ax1.imshow(data,origin='bottom',cmap=cm.gray,interpolation='none')
	plt.show()
