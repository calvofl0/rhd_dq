#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
This nMBP module locates non-magnetic bright points using an indicator based
on three parameters: temperature, vorticity and depth at isosurface tau=1.

This indicator is a growing function of its parameters. Local maxima above
some threshold of this indicator are considered to be nMBPs. The indicator
has been experimentally built, and hasn't any theoretical justification.

When nMBPs are too close, only the brightest one (in the intensity map) is
retained, considering that the underlying structure is that of one single
nMBP providing two close local maxima in the indicator.

Also if the position of the hypothetical nMBP in the intensity map is too
far away from the underlying structure, the nMBP is eliminated.
'''

import numpy as np
from slicingTools import vorticity
from snake import find_min, select_disk, d2
from analyze2d import getGranules, circleFootprint

import pickle
with open('rad1145.pickle') as f:
	rad=pickle.load(f)
with open('tau1145.pickle') as f:
	tau1=pickle.load(f)
with open('snap1145.pickle') as f:
	v1,v2,v3,rho=pickle.load(f)
	v=vorticity(v1,v2)
with open('vorticity1145.pickle') as f:
	vtau1=pickle.load(f)
with open('T1145.pickle') as f:
	T=pickle.load(f)

normalize = lambda data: (data-np.min(data))/(np.max(data)-np.min(data))

def computeIndicator(T, tau1, vorticity):
	return np.sqrt(normalize(T)**2*normalize(-tau1)**2*normalize(np.abs(vorticity)**2))

def nearestMin(data, p, max_dist=20):
	mins = find_min(data)
	dist = []
	intensity = []
	for m in mins:
		dist.append(d2(p,m))
		intensity.append(data[int(m[0]),int(m[1])])
	sortTable = np.argsort(np.array(intensity))
	intensity = np.array(intensity)[sortTable]
	mins = np.array(mins)[sortTable]
	dist = np.array(dist)[sortTable]
	for m,d in zip(mins,dist):
		if d <= max_dist**2: return (int(m[0]), int(m[1])), 0
	# print('Warning: nearest min is further than max_dist.')
	m = mins[np.argmin(dist)]
	return (int(m[0]), int(m[1])), 1

def circleFootprint(N):
	fp = np.zeros((N,N), dtype=bool)
	r = (N-1.)/2.
	for i in range(N):
		for j in range(N):
			if (i-r)**2+(j-r)**2 <= r**2: fp[i,j] = True
	return fp

def getNMBPs(indicator, intensity, granules=None, l=50, footprint=None):
	if granules == None: granules = getGranules(intensity, shrink=0.)
	if footprint == None: fp = circleFootprint(l)
	elif type(footprint) == int: fp = circleFootprint(footprint)
	else: fp = footprint
	p = find_min(-indicator)
	p = [(int(p0[0]), int(p0[1])) for p0 in p]
	p = [p0 for p0 in p if indicator[p0]>.2]
	zones = []
	toRemove = set()
	#minimas = []
	nmbps = []
	for p00,p01 in p:
		nx, ny = np.shape(indicator)
		mgrid=np.meshgrid(np.arange(p00-l,p00+l+1)%nx,np.arange(p01-l,p01+l+1)%ny)
		mgrid=(mgrid[0].T,mgrid[1].T)
		q,err = nearestMin(-intensity[mgrid],(l+1,l+1))
		if err != 0:
			print('Warning: nearest nMBP in intensity map to point ('+str(p00)+', '+str(p01)+') is further than max_dist.')
		s,d,b = select_disk(-intensity[mgrid],q)
		zones.append(zip((s[0]+p00-l)%nx,(s[1]+p01-l)%ny))
		nmbps.append(((q[0]+p00-l)%nx,(q[1]+p01-l)%ny))
	for i in range(len(p)):
		if p[i] not in zones[i]:
			toRemove.add(i)
			continue
		for j in range(len(p)):
			if j>=i: break
			s = set.intersection(set(zones[i]), set(zones[j]))
			if len(s) == 0: continue
			if len(zones[j])>len(zones[i]): toRemove.add(i)
			else: toRemove.add(j)
	if len(toRemove)>0:
		toRemove = np.sort(np.array(list(toRemove)))-np.arange(len(toRemove))
	for i in toRemove:
		p.pop(i)
		zones.pop(i)
		nmbps.pop(i)
	background = []
	contrast_local = []
	contrast_global = []
	int_global = np.mean(intensity)
	for i in range(len(nmbps)):
		xshift = int(nmbps[i][0]-(np.size(fp, axis=0)-1)/2.)
		yshift = int(nmbps[i][1]-(np.size(fp, axis=1)-1)/2.)
		NX, NY = np.shape(intensity)
		n = np.count_nonzero(fp)
		fp_mask = np.array(np.where(fp))+np.repeat([[xshift],[yshift]], n, axis=1)
		fp_mask = (np.mod(fp_mask[0], NX), np.mod(fp_mask[1], NY))
		nmbp_mask = zip(*zones[i])
		mask = zip(*set.difference(set(zip(*fp_mask)), set(zip(*nmbp_mask))))
		dnh = intensity[mask]*(1-granules[mask]) # Dark neighbourhood
		background.append(np.sum(dnh)/np.count_nonzero(dnh))
		contrast_local.append(np.mean(intensity[zip(*zones[i])]) / background[i])
		contrast_global.append(np.mean(intensity[zip(*zones[i])]) / int_global)
	
	return zip(p, np.array(contrast_local)-1., np.array(contrast_global)-1.)
