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

Via the functions export_nMBPs, import_nMBPs, importNext, report and
setActive the nMBPs can be studies. The typical strategy is:

1. export_nMBPs(fullfile, parfile, meanfile, height)
	-> This will create uio files with cuts centred on each nMBP and
	   at height 'height'. It also creates an index file '.nmbps'
	   This step does not requires matplotlib

2. import_nMBPs(indexfile, parfile)
	-> This loads the header of the index file, and these next steps
	   should be done in a computer where matplotlib is installed

3. importNext()
	-> This should be called as many times as there are nMBPs. At each
	   call a new nMBP is loaded. It returns a list with:
		* The coordinates (a tuple)
		* Contrast (with respect to local neighbourhood)
		* Contrast (with respect to global average)
		* Diameter of the nMBP
		* A flag specifying wheter or not the nMBP is set active
		  (this flag is switched off for non-valid nMBPs after
		   visual examination)

4. report()
	-> Prints a report of the physical properties of the current nMBP
	   together with the corresponding plots

5. setActive(active_flag)
	-> Set the active_flag of the current nMBP (this modifies the index)

6. close()
	-> Closes the index

'''

import numpy as np
import pybold
import snake
import xdrlib
from slicingTools import vorticity
from snake import find_min, select_disk, d2
from analyze2d import getGranules, circleFootprint
from analyze3d import wilson_depression, contrast
from pybold import uio_struct, varAtLevel
try:
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	from matplotlib.colors import colorConverter
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	from analyze3d import plotv_slice
except ImportError: print('Warning: matplotlib is missing.') 

if not '_model' in globals(): _model = uio_struct()
if not '_modelfile' in globals(): _modelfile = None
if not '_parfile' in globals(): _parfile = None
if not '_meanfile' in globals(): _meanfile = None
if not '_indexfile' in globals(): _indexfile = None
if not '_index' in globals(): _index = None

if not '_tau1' in globals(): _tau1 = None
if not '_vtau1' in globals(): _vtau1 = None
if not '_T' in globals(): _T = None

if not '_allrad' in globals(): _allrad = None
if not '_rad' in globals(): _rad = None

if not '_xdr' in globals(): _xdr = None
if not '_Nx' in globals(): _Nx = 0
if not '_Ny' in globals(): _Ny = 0
if not '_Npts' in globals(): _Npts = 0
if not '_pt' in globals(): _pt = 0
if not '_parameters' in globals(): _parameters = None
if not '_pos' in globals(): _pos = 0

def load_uio(modelfile, parfile=None, meanfile=None):
	'''
	Loads an uio file.
	'''
	global _model, _modelfile, _parfile, _meanfile, _allrad, _rad
	#_model = uio_struct()
	_modelfile = modelfile
	_parfile   = None
	if meanfile != None:
		_meanmodel = uio_struct()
		_meanmodel.open(meanfile)
		_meanfile = meanfile
		_meanmodel.header
		_allrad = dict()
		while _meanmodel.next >= 0:
			_allrad[_meanmodel.modelitime] = np.copy(_meanmodel.rad.intb3_r[:,:,0])
		_meanmodel.close
	if parfile:
		_model.open(modelfile, parfile)
		_parfile = parfile
	else: _model.open(modelfile)
	_model.header
	return next()

def next():
	'''
	Loads next snapshot
	'''
	global _model, _modelfile, _parfile, _meanfile, _allrad, _rad
	if _modelfile is None: return -1
	print('Get next snapshot')
	if _model.next < 0: return -1
	print('Got next snapshot')
	if _meanfile != None:
		itime = _model.modelitime
		if _allrad.has_key(itime):
			_rad = _allrad[itime]
		else: print('Warning: snapshot not found in MEAN file.')
	return _model.modelitime

def close():
	'''
	Closes model or index file depending on the context
	'''
	global _model, _index
	if _index != None:
		_index.close()
		_index = None
	else:
		_model.close

def get_tau1vars(model=_model):
	'''
	Compute physical quantities that will be used as indicators
	'''
	global _tau1, _vtau1, _T
	# tau1
	NX, NY, NZ = np.size(model.z.xc1), np.size(model.z.xc2), np.size(model.z.xc3)
	Z = np.repeat(np.array([np.repeat([model.z.xc3[0,0,:]],NY,axis=0)]),NX,axis=0)
	_tau1 = varAtLevel(Z,model.dq.tau,1.)
	# vtau1
	v1, v2 = model.z.v1, model.z.v2
	#dxb1 = model.z.xb1[1:,:,:] - model.z.xb1[:-1,:,:]
	#dxb2 = model.z.xb2[:,1:,:] - model.z.xb2[:,:-1,:]
	#dxb1 = np.repeat(dxb1, NY, 1)
	#dxb1 = np.repeat(dxb1, NZ, 2)
	#dxb2 = np.repeat(dxb2, NX, 0)
	#dxb2 = np.repeat(dxb2, NZ, 2)
	#dv1d2 = (np.roll(v1, -1, 1) - v1) / dxb2
	#dv2d1 = (np.roll(v2, -1, 0) - v2) / dxb1
	#_vtau1 = varAtLevel(dv1d2-dv2d1,model.dq.tau,1.)
	v1 = varAtLevel(v1, model.dq.tau, 1.)
	v2 = varAtLevel(v2, model.dq.tau, 1.)
	_vtau1 = vorticity(v1, v2)
	# T
	_T = varAtLevel(model.dq.T,model.dq.tau,1.)

#import pickle
#with open('rad1145.pickle') as f:
#	rad=pickle.load(f)
#with open('tau1145.pickle') as f:
#	tau1=pickle.load(f)
#with open('snap1145.pickle') as f:
#	v1,v2,v3,rho=pickle.load(f)
#	v=vorticity(v1,v2)
#with open('vorticity1145.pickle') as f:
#	vtau1=pickle.load(f)
#with open('T1145.pickle') as f:
#	T=pickle.load(f)

normalize = lambda data: (data-np.min(data))/(np.max(data)-np.min(data))

def computeIndicator(T=None, tau1=None, vorticity=None):
	'''
	This functions computes the indicator used to select nMBPs
	It is the core of the detection algorithm
	'''
	if T is None: T=_T
	if tau1 is None: tau1=_tau1
	if vorticity is None: vorticity=_vtau1
	return np.sqrt(normalize(T)**2*normalize(-tau1)**2*normalize(np.abs(vorticity)**2))

def nearestMin(data, p, max_dist=20):
	'''
	Returns the nearest minimum to p together to a flag checking whether
	or not it is not further away than max_dist
	'''
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
	'''
	Creates a circular "footprint" (or stencil) of diameter N
	'''
	fp = np.zeros((N,N), dtype=bool)
	r = (N-1.)/2.
	for i in range(N):
		for j in range(N):
			if (i-r)**2+(j-r)**2 <= r**2: fp[i,j] = True
	return fp

def getNMBPs(indicator, intensity=None, granules=None, l=50, footprint=None):
	'''
	Gets all nMBPs using the given indicator. If an intensity map is not
	provided, temperature at isosurface tau=1. is taken in place. If
	granules were already computed, they can be provided here to avoid
	a new computations. If footprint is not provided a circular footprint
	(or stencil) of diameter l is used.
	'''
	global _model, _T
	if intensity is None: intensity = _T
	if granules is None and _model is None: granules = getGranules(intensity, shrink=0.)
	if granules is None and _model != None:
		v3 = varAtLevel(_model.z.v3, _model.dq.tau, 1.)
		granules = v3>0.
	if footprint is None: fp = circleFootprint(l)
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
	diametre = []
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
		diametre.append(2.*np.sqrt(len(zones[i])/np.pi))
	
	return zip(p, np.array(contrast_local)-1., np.array(contrast_global)-1., diametre)


def export_nMBPs(modelfile, parfile, meanfile, height, box=(100,100,100)):
	'''
	Export all nMBPs to smaller boxes of size 'box' centred at the nMBPs
	position and height 'height' and saves the boxes to uio files.
	An index file is further written, including intensity maps.
	'''
	global _model, _modelfile, _parfile, _meanfile, _allrad, _rad
	itime = load_uio(modelfile, parfile, meanfile)

	parts = modelfile.split('.')
	l = max(0, len(parts)-1)
	filename = ''.join(parts[:l])

	xdr = xdrlib.Packer()
	xdr.pack_string('nmbps_index')
	Nx = _model.z.dimension[1][0]-_model.z.dimension[0][0]+1
	Ny = _model.z.dimension[1][1]-_model.z.dimension[0][1]+1
	xdr.pack_int(Nx)
	xdr.pack_int(Ny)
	
	while itime >= 0:
		xdr.pack_bool(True)
		xdr.pack_array(_rad.flatten(), xdr.pack_float)
		get_tau1vars()
		indicator = computeIndicator()
		nmbps = getNMBPs(indicator, intensity=_rad)
		xdr.pack_int(len(nmbps))
		for pt in nmbps:
			xdr.pack_int(_model.modelitime)
			pX = pt[0][0]
			pY = pt[0][1]
			xdr.pack_int(pt[0][0])
			xdr.pack_int(pt[0][1])
			xdr.pack_float(pt[1])
			xdr.pack_float(pt[2])
			xdr.pack_float(pt[3])
			xdr.pack_bool(True)
			# Save uio file here!
			_model.export(filename+'_'+str(_model.modelitime)+'_'+str(pX)+'_'+str(pY)+'.uio', centre=(pX, pY, height), box=box)
		itime = next()
	xdr.pack_bool(False)
	with open(filename+'.nmbps', 'wb') as f:
		f.write(xdr.get_buffer())

def show_nMBPs(indexfile):
	'''
	Show all available nMBPs in indexfile
	'''

	with open(indexfile, 'rb') as f:
		xdr = xdrlib.Unpacker(f.read())
	if xdr.unpack_string() != 'nmbps_index':
		print('Error: unknown file.')
		return
	Nx = xdr.unpack_int()
	Ny = xdr.unpack_int()
	pt = 0
	Npts = 0
	while xdr.unpack_bool():
		rad = np.array(xdr.unpack_array(xdr.unpack_float)).reshape((Nx, Ny))
		Npts = xdr.unpack_int()
		pt = 0
		pX = []
		pY = []
		circs = []
		while pt < Npts:
			pt += 1
			modelitime = xdr.unpack_int()
			pX = xdr.unpack_int()
			pY = xdr.unpack_int()
			contrast_local  = xdr.unpack_float()
			contrast_global = xdr.unpack_float()
			diametre        = xdr.unpack_float()
			active_flag	= xdr.unpack_bool()
			if active_flag:
				circs += [plt.Circle((pX,pY),50,color='r',fill=False)]
		fig = plt.figure()
		plt.imshow(rad.T,origin='bottom',cmap=plt.cm.gray)
		for c in circs:
			fig.gca().add_artist(c)

def import_nMBPs(indexfile, parfile):
	'''
	Loads nMBPs that where written by the corresponding exporting
	function. Actually this only reads the header of the index file.
	'''
	global _model, _modelfile, _parfile, _indexfile, _xdr, _Nx, _Ny, _Npts, _pt
	_parfile = parfile
	_indexfile = indexfile.strip('.nmbps')
	with open(indexfile, 'rb') as f:
		_xdr = xdrlib.Unpacker(f.read())
	if _xdr.unpack_string() != 'nmbps_index':
		print('Error: unknown file.')
		return
	_Nx = _xdr.unpack_int()
	_Ny = _xdr.unpack_int()
	_pt = 0
	_Npts = 0

def importNext():
	'''
	Loads the next available nMBP (after calling import_nMBPs)
	'''
	global _model, _modelfile, _parfile, _indexfile, _rad, _xdr, _Nx, _Ny, _Npts, _pt, _parameters, _pos
	if _pt < 0:
		print('No more nMBPs available.')
		return
	if _pt >= _Npts-1:
		if not _xdr.unpack_bool():
			_pt = -1
			print('No more nMBPs available.')
			return
		_rad = np.array(_xdr.unpack_array(_xdr.unpack_float)).reshape((_Nx, _Ny))
		_Npts = _xdr.unpack_int()
		_pt = 0
		if _Npts <= _pt:
			_pt = -1
			print('No more nMBPs available.')
			return
	else:
		_pt = _pt+1

	modelitime = _xdr.unpack_int()
	pX = _xdr.unpack_int()
	pY = _xdr.unpack_int()
	contrast_local  = _xdr.unpack_float()
	contrast_global = _xdr.unpack_float()
	#dpos		= _xdr.get_position()
	diametre        = _xdr.unpack_float()
	_pos		= _xdr.get_position()
	active_flag	= _xdr.unpack_bool()
	#with open(_indexfile+'.nmbps', 'r+b') as f:
	#	f.seek(dpos)
	#	xdr = xdrlib.Packer()
	#	xdr.pack_float(2.*diametre)
	#	f.write(xdr.get_buffer())
	#return True
	# Load uio_file here!
	_modelfile = _indexfile + '_' + str(modelitime) + '_' + str(pX) + '_' + str(pY) + '.uio'
	if load_uio(_modelfile, _parfile) >= 0:
		_model.close

	_parameters = [(pX, pY), contrast_local, contrast_global, diametre, active_flag]
	return _parameters

def setActive(active_flag):
	'''
	Set the active_flag of the current nMBP. This modified the index file.
	'''
	global _indexfile, _pos
	_parameters[4] = active_flag
	with open(_indexfile+'.nmbps', 'r+b') as f:
		f.seek(_pos)
		xdr = xdrlib.Packer()
		xdr.pack_bool(active_flag)
		f.write(xdr.get_buffer())

def report():
	'''
	Prints a report of the current nMBP together with the corresponding
	plots.
	'''
	global _model, _parameters, _rad
	model = _model
	tau = 1.
	r = 3
	tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
	tau1=int(round(np.mean(tau1_level)))
	tau1min=int(round(np.min(tau1_level)))
	seed=[_parameters[0][0]-model.z.dimension[0][0],_parameters[0][1]-model.z.dimension[0][1],tau1]
	s=snake.snake_from_box(model.z.rho,radius=r,start=tau1,periodic=False,seed=seed)
	tau1_sindex = int(list(s[:,2]).index(tau1))
	xmean = int(s[tau1_sindex,0])
	ymean = int(s[tau1_sindex,1])
	odb=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean),threshold=.5)
	rot_index = np.argmin(model.z.v1[:,:,tau1][odb[0]]**2)
	xrot,yrot = odb[0][0][rot_index], odb[0][1][rot_index]
	xc1, xc2 = model.z.xc1.flatten()/1.e5, model.z.xc2.flatten()/1.e5
	xb1, xb2 = model.z.xb1.flatten()/1.e5, model.z.xb2.flatten()/1.e5
	ext = [xc1[0],xc1[-1],xc2[0],xc2[-1]]
	#Ttau1 = pybold.varAtLevel(model.dq.T,model.dq.tau,1.)
	rhotau1 = pybold.varAtLevel(model.z.rho,model.dq.tau,1.)
	#plt.imshow(model.z.rho[:,:,tau1].T,origin='bottom',cmap=plt.cm.gray,extent=ext)
        plt.imshow(rhotau1.T,origin='bottom',extent=ext)
	plt.xlim(xc1[0],xc1[-1])
	plt.ylim(xc2[0],xc2[-1])
	plt.xlabel('Spatial horizontal position (X axis) [km]')
	plt.ylabel('Spatial horizontal position (Y axis) [km]')
        oring=np.zeros_like(model.z.rho[:,:,tau1],dtype=bool)
        oring[odb[0]]=True
        plt.imshow(oring.T,origin='bottom',alpha=.15,interpolation='none',extent=ext,cmap=plt.cm.gray_r)
	X, Y = np.meshgrid(xc1, xc2)
	arrfreq = 5
	v1 = pybold.varAtLevel(model.z.v1, model.dq.tau, 1.)
	v2 = pybold.varAtLevel(model.z.v2, model.dq.tau, 1.)
	#plt.quiver(X[::arrfreq,::arrfreq], Y[::arrfreq,::arrfreq], model.z.v1[::arrfreq,::arrfreq,tau1].T,model.z.v2[::arrfreq,::arrfreq,tau1].T)
	plt.quiver(X[::arrfreq,::arrfreq], Y[::arrfreq,::arrfreq], v1[::arrfreq,::arrfreq].T,v2[::arrfreq,::arrfreq].T)
	fmt='{0: <20}'
        if _modelfile: print(fmt.format('File:')+_modelfile)
        else: print(fmt.format('File:')+'Unknown')
	wd=wilson_depression(model,s,odb,tau,r)
        print(fmt.format('Wilson depression:')+str(wd))
        c_rho=contrast(model.z.rho,model,s,odb,tau,r)
        print(fmt.format('Density contrast:')+str(c_rho))
        c_p=contrast(model.dq.P,model,s,odb,tau,r)
        print(fmt.format('Pressure contrast:')+str(c_p))
	if _modelfile:
                title = _modelfile.split('.')
                title = '.'.join(title[:max(1,len(title)-1)])
	boxtext = 'Wilson depression: '+str(int(round(wd/1.e5)))+' [km]\n'
        boxtext+= 'Density contrast: '+str(int(round(100*c_rho[1])))+'%, '
        boxtext+= str(int(round(100*c_rho[0])))+'%\n'
        boxtext+= 'Pressure contrast: '+str(int(round(100*c_p[1])))+'%, '
        boxtext+= str(int(round(100*c_p[0])))+'%'
	diam = _parameters[3]
	delta = model.z.xc1[1,:,:]-model.z.xc1[0,:,:]
        boxtext+= '\nDiametre: '+str(int(round(diam*delta/1.e5)))+' [km]'
	c_I = (_parameters[1], _parameters[2])
        boxtext+= '\nIntensity contrast: '+str(int(round(100*c_I[0])))+'%, '
        boxtext+= str(int(round(100*c_I[1])))+'%'
	boxtext+= '\nLifetime: '+str(int(round(120.)))+' [s]'
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	deltaX = xb1[1]-xb1[0]
	deltaY = xb2[1]-xb2[0]
	xlim_min = xc1[0]-model.z.dimension[0][0]*deltaX
	ylim_min = xc2[0]-model.z.dimension[0][1]*deltaY
	szX, szY = np.shape(_rad.T)
	ext = [xlim_min,xlim_min+szX*deltaX,ylim_min,ylim_min+szY*deltaY]
	ax.imshow(_rad.T, origin='bottom', cmap=plt.cm.gray, extent=ext)
	ax.set_xlim(xlim_min,xlim_min+szX*deltaX)
	ax.set_ylim(ylim_min,ylim_min+szY*deltaY)
	ax.set_xlabel('Spatial horizontal position (X axis) [km]')
	ax.set_ylabel('Spatial horizontal position (Y axis) [km]')
	c1_idx = _parameters[0][0] - model.z.dimension[0][0]
	c2_idx = _parameters[0][1] - model.z.dimension[0][1]
	centre = [xc1[c1_idx], xc2[c2_idx]]
	c_r    = 50.*np.sqrt(deltaX**2+deltaY**2)/np.sqrt(2.)
	circ = plt.Circle(centre, radius=c_r, color='r', fill=False)
	ax.add_patch(circ)
	sY,sX=plotv_slice(model.z.rho,model,s,dq=True,tau=tau,r=r,boxtext=boxtext,show=True,tight=False)
        plt.close()
