#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a high-level Python module to analyse individual non-magnetic bright
points. It uses the lower-level "snake" module and the pybold module for
I/O.

Non-magnetic bright points should be stored individually in small uio files.
This module provides the following functions (optional arguments are marked
with '*' and are derived automatically if not specified and if the model is
loaded with the load_uio function rather than directly with pybold):

	- nmbp = load_uio(modelfile, parfile*): loads an uio file in nmbp
	- grad = gradient(arr, nmbp*): computes gradient of arr from a nMBP.
	- a_c  = a_centripetal(nmbp*,s*,r*): computes centripetal
	         acceleration from rotation around column s (s is given by
		 the "snake" module)
	- a_g  = a_gradient(nmbp*): computes the actual acceleration due to
	         gradient pressure
	- wd   = wilson_depression(model*, s*, odb*, tau*, r*): computes the
	         Wilson depression
	- l,g  = contrast(arr, model*, s*, odb*, tau*, r*): computes local
	         and global contrast of arr; (arr_i-<arr>)/<arr>. <arr> is
		 taken in the local neighbourhood and on the whole z=cst
		 slice, at <tau>=tau (=1 by default).
	- Y,X  = plotv_slice(arr, model*, s*, dq*, tau*, r*): plots vertical
		 slices of arr at y=cst and x=cst respectively. Depending on
		 optional parameters, also plots tau=tau (=1 by default)
		 isosurface. Marks in red the "eye" of the nMBP.
"""

import numpy as np
try:
	from scipy.interpolate import griddata, interp1d
except ImportError: print('Warning: scipy is missing.')
try:
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	#import matplotlib.colors as colors
except ImportError: print('Warning: matplotlib is missing.')
import pybold
import snake
from savitzky_golay import savitzky_golay

if not '_model' in globals(): _model = pybold.uio_struct()
if not '_modelfile' in globals(): _modelfile = None
if not '_parfile' in globals(): _parfile = None

def load_uio(modelfile, parfile=None):
	'''
	Loads an uio file.
	'''
	global _model, _modelfile, _parfile
	_model = pybold.uio_struct()
	_modelfile = modelfile
	_parfile = None
	if parfile:
		_model.open(modelfile, parfile)
		_parfile = parfile
	else: _model.open(modelfile)
	_model.header
	_model.next
	return _model

def meshgrid(model=None):
	'''
	3D "meshgrid", returns X, Y, Z so that
		* X_ijk = i
		* Y_ijk = j
		* Z_ijk = k
	'''
	if not model: model=_model
	NX, NY, NZ = np.size(model.z.xc1), np.size(model.z.xc2), np.size(model.z.xc3)
	X = np.repeat(np.array([np.repeat([model.z.xc1[:,0,0]],NY,axis=0)]),NZ,axis=0).T
	Y = np.repeat(np.array([np.repeat([model.z.xc2[0,:,0]],NZ,axis=0).T]),NX,axis=0)
	Z = np.repeat(np.array([np.repeat([model.z.xc3[0,0,:]],NY,axis=0)]),NX,axis=0)
	return X,Y,Z

def meshgrid2D(model=None):
	'''
	2D "meshgrid", returns X, Y so that
		* X_ij = i
		* Y_ij = j
	'''
	if not model: model=_model
	X, Y = np.meshgrid(model.z.xc1[:,0,0],model.z.xc2[0,:,0])
	return X.T, Y.T

def gradient(arr, model=None, z=None):
	'''
	Computes gradient of arr. If z is give, returns a 2d array along
	each dimension.	If z is not given, returns arrays of the same size
	than arr for each dimension of arr.
	'''
	if not model: model=_model
	s = np.shape(arr)
	dtp = arr.dtype
	X, Y, Z = meshgrid(model)
	if z != None:
		gradX = np.empty((s[0],s[1]), dtype=dtp)
		gradY = np.empty((s[0],s[1]), dtype=dtp)
		gradZ = np.empty((s[0],s[1]), dtype=dtp)
		gradX[1:-1,:] = (arr[2:,:,z]-arr[:-2,:,z])/(X[2:,:,z]-X[:-2,:,z])
		gradY[:,1:-1] = (arr[:,2:,z]-arr[:,:-2,z])/(Y[:,2:,z]-Y[:,:-2,z])
		if z>0 and z<s[2]-1:
			gradZ[:,:] = (arr[:,:,z+1]-arr[:,:,z-1])/(Z[:,:,z+1]-Z[:,:,z-1])
		elif z==0:
			gradZ[:,:] = arr[:,:,1]+(arr[:,:,1]-arr[:,:,2])/(Z[:,:,1]-Z[:,:,2])*(Z[:,:,0]-Z[:,:,1])
		elif z==s[2]-1:
			gradZ[:,:] = arr[:,:,s[2]-1]+(arr[:,:,s[2]-1]-arr[:,:,s[2]-2])/(Z[:,:,s[2]-1]-Z[:,:,s[2]-2])*(Z[:,:,s[2]]-Z[:,:,s[2]-1])
		gradX[0,:] = gradX[1,:]+(gradX[1,:]-gradX[2,:])/(X[1,:,z]-X[2,:,z])*(X[0,:,z]-X[1,:,z])
		gradX[-1,:] = gradX[-2,:]+(gradX[-2,:]-gradX[-3,:])/(X[-2,:,z]-X[-3,:,z])*(X[-1,:,z]-X[-2,:,z])
		gradY[:,0] = gradY[:,1]+(gradY[:,1]-gradY[:,2])/(Y[:,1,z]-Y[:,2,z])*(Y[:,0,z]-Y[:,1,z])
		gradY[:,-1] = gradY[:,-2]+(gradY[:,-2]-gradY[:,-3])/(Y[:,-2,z]-Y[:,-3,z])*(Y[:,-1,z]-Y[:,-2,z])
		return gradX, gradY, gradZ
	else:
		gradX = np.empty_like(arr)
		gradY = np.empty_like(arr)
		if len(s) == 3:
			gradZ = np.empty_like(arr)
			gradX[1:-1,:,:] = (arr[2:,:,:]-arr[:-2,:,:])/(X[2:,:,:]-X[:-2,:,:])
			gradY[:,1:-1,:] = (arr[:,2:,:]-arr[:,:-2,:])/(Y[:,2:,:]-Y[:,:-2,:])
			gradZ[:,:,1:-1] = (arr[:,:,2:]-arr[:,:,:-2])/(Z[:,:,2:]-Z[:,:,:-2])
			gradX[0,:,:] = gradX[1,:,:]+(gradX[1,:,:]-gradX[2,:,:])/(X[1,:,:]-X[2,:,:])*(X[0,:,:]-X[1,:,:])
			gradX[-1,:,:] = gradX[-2,:,:]+(gradX[-2,:,:]-gradX[-3,:,:])/(X[-2,:,:]-X[-3,:,:])*(X[-1,:,:]-X[-2,:,:])
			gradY[:,0,:] = gradY[:,1,:]+(gradY[:,1,:]-gradY[:,2,:])/(Y[:,1,:]-Y[:,2,:])*(Y[:,0,:]-Y[:,1,:])
			gradY[:,-1,:] = gradY[:,-2,:]+(gradY[:,-2,:]-gradY[:,-3,:])/(Y[:,-2,:]-Y[:,-3,:])*(Y[:,-1,:]-Y[:,-2,:])
			gradZ[:,:,0] = gradZ[:,:,1]+(gradZ[:,:,1]-gradZ[:,:,2])/(Z[:,:,1]-Z[:,:,2])*(Z[:,:,0]-Z[:,:,1])
			gradZ[:,:,-1] = gradZ[:,:,-2]+(gradZ[:,:,-2]-gradZ[:,:,-3])/(Z[:,:,-2]-Z[:,:,-3])*(Z[:,:,-1]-Z[:,:,-2])
			return gradX, gradY, gradZ
		else:
			gradX[1:-1,:] = (arr[2:,:]-arr[:-2,:])/(X[2:,:]-X[:-2,:])
			gradY[:,1:-1] = (arr[:,2:]-arr[:,:-2])/(Y[:,2:]-Y[:,:-2])
			gradX[0,:] = gradX[1,:]+(gradX[1,:]-gradX[2,:])/(X[1,:]-X[2,:])*(X[0,:]-X[1,:])
			gradX[-1,:] = gradX[-2,:]+(gradX[-2,:]-gradX[-3,:])/(X[-2,:]-X[-3,:])*(X[-1,:]-X[-2,:])
			gradY[:,0] = gradY[:,1]+(gradY[:,1]-gradY[:,2])/(Y[:,1]-Y[:,2])*(Y[:,0]-Y[:,1])
			gradY[:,-1] = gradY[:,-2]+(gradY[:,-2]-gradY[:,-3])/(Y[:,-2]-Y[:,-3])*(Y[:,-1]-Y[:,-2])
			return gradX, gradY

def gradient2(arr, model=None, dq=None, tau=None, z=None):
	'''
	Computes gradient of arr from a nMBP.
	'''
	if not model: model=_model
	dim = len(np.shape(arr))
	if dim == 3 and tau != None: to2d=True
	else: to2d=False
	if tau == None: tau=1.
	if dq == None: dq=bool(_parfile)
	if dq:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
	else: tau1=int(np.size(model.z.xc3)/2.)
	if z != None: tau1=z
	if to2d:
		gradX = np.empty_like(arr[:,:,tau1])
		gradY = np.empty_like(arr[:,:,tau1])
		if dim == 3:
			gradZ = np.empty_like(arr[:,:,tau1])
	else:
		gradX = np.empty_like(arr)
		gradY = np.empty_like(arr)
		if dim == 3:
			gradZ = np.empty_like(arr)
	[[m1,m2,m3],[n1,n2,n3]] = model.z.dimension
	X, Y, Z = meshgrid(model)
	#if dim == 3: X, Y, Z = meshgrid(model)
	#else: X, Y = meshgrid2D(model)
	if to2d:
		gradX[1:-1,:] = (arr[2:,:,tau1]-arr[:-2,:,tau1])/(X[2:,:,tau1]-X[:-2,:,tau1])
		gradY[:,1:-1] = (arr[:,2:,tau1]-arr[:,:-2,tau1])/(Y[:,2:,tau1]-Y[:,:-2,tau1])
		if tau1>0 and tau1<n3-m3:
			gradZ[:,:] = (arr[:,:,tau1+1]-arr[:,:,tau1-1])/(Z[:,:,tau1+1]-Z[:,:,tau1-1])
		elif tau1==0:
			gradZ[:,:] = arr[:,:,1]+(arr[:,:,1]-arr[:,:,2])/(Z[:,:,1]-Z[:,:,2])*(Z[:,:,0]-Z[:,:,1])
		elif tau1==n3-m3:
			gradZ[:,:] = arr[:,:,n3-m3-1]+(arr[:,:,n3-m3-1]-arr[:,:,n3-m3-2])/(Z[:,:,n3-m3-1]-Z[:,:,n3-m3-2])*(Z[:,:,n3-m3]-Z[:,:,n3-m3-1])
		gradX[0,:] = arr[1,:,tau1]+(arr[1,:,tau1]-arr[2,:,tau1])/(X[1,:,tau1]-X[2,:,tau1])*(X[0,:,tau1]-X[1,:,tau1])
		gradX[-1,:] = arr[-2,:,tau1]+(arr[-2,:,tau1]-arr[-3,:,tau1])/(X[-2,:,tau1]-X[-3,:,tau1])*(X[-1,:,tau1]-X[-2,:,tau1])
		gradY[:,0] = arr[:,1,tau1]+(arr[:,1,tau1]-arr[:,2,tau1])/(Y[:,1,tau1]-Y[:,2,tau1])*(Y[:,0,tau1]-Y[:,1,tau1])
		gradY[:,-1] = arr[:,-2,tau1]+(arr[:,-2,tau1]-arr[:,-3,tau1])/(Y[:,-2,tau1]-Y[:,-3,tau1])*(Y[:,-1,tau1]-Y[:,-2,tau1])
		return (gradX, gradY, gradZ)
	elif dim == 3:
		gradX[1:-1,:,:] = (arr[2:,:,:]-arr[:-2,:,:])/(X[2:,:,:]-X[:-2,:,:])
		gradY[:,1:-1,:] = (arr[:,2:,:]-arr[:,:-2,:])/(Y[:,2:,:]-Y[:,:-2,:])
		gradZ[:,:,1:-1] = (arr[:,:,2:]-arr[:,:,:-2])/(Z[:,:,2:]-Z[:,:,:-2])
		gradX[0,:,:] = gradX[1,:,:]+(gradX[1,:,:]-gradX[2,:,:])/(X[1,:,:]-X[2,:,:])*(X[0,:,:]-X[1,:,:])
		gradX[-1,:,:] = gradX[-2,:,:]+(gradX[-2,:,:]-gradX[-3,:,:])/(X[-2,:,:]-X[-3,:,:])*(X[-1,:,:]-X[-2,:,:])
		gradY[:,0,:] = gradY[:,1,:]+(gradY[:,1,:]-gradY[:,2,:])/(Y[:,1,:]-Y[:,2,:])*(Y[:,0,:]-Y[:,1,:])
		gradY[:,-1,:] = gradY[:,-2,:]+(gradY[:,-2,:]-gradY[:,-3,:])/(Y[:,-2,:]-Y[:,-3,:])*(Y[:,-1,:]-Y[:,-2,:])
		gradZ[:,:,0] = gradZ[:,:,1]+(gradZ[:,:,1]-gradZ[:,:,2])/(Z[:,:,1]-Z[:,:,2])*(Z[:,:,0]-Z[:,:,1])
		gradZ[:,:,-1] = gradZ[:,:,-2]+(gradZ[:,:,-2]-gradZ[:,:,-3])/(Z[:,:,-2]-Z[:,:,-3])*(Z[:,:,-1]-Z[:,:,-2])
		return (gradX, gradY, gradZ)
	else:
		gradX[1:-1,:] = (arr[2:,:]-arr[:-2,:])/(X[2:,:]-X[:-2,:])
		gradY[:,1:-1] = (arr[:,2:]-arr[:,:-2])/(Y[:,2:]-Y[:,:-2])
		gradX[0,:] = gradX[1,:]+(gradX[1,:]-gradX[2,:])/(X[1,:]-X[2,:])*(X[0,:]-X[1,:])
		gradX[-1,:] = gradX[-2,:]+(gradX[-2,:]-gradX[-3,:])/(X[-2,:]-X[-3,:])*(X[-1,:]-X[-2,:])
		gradY[:,0] = gradY[:,1]+(gradY[:,1]-gradY[:,2])/(Y[:,1]-Y[:,2])*(Y[:,0]-Y[:,1])
		gradY[:,-1] = gradY[:,-2]+(gradY[:,-2]-gradY[:,-3])/(Y[:,-2]-Y[:,-3])*(Y[:,-1]-Y[:,-2])
		return (gradX, gradY)

def advection(arr, model=None, z=None):
	'''
	Computes advection of arr.
	'''
	if not model: model=_model
	if z != None:
		v1,v2,v3 = model.z.v1[:,:,z],model.z.v2[:,:,z],model.z.v3[:,:,z]
		#v1,v2,v3 = b1[:,:,z],b2[:,:,z],b3[:,:,z]
	else: v1,v2,v3 = model.z.v1,model.z.v2,model.z.v3
	#	v1,v2,v3 = b1,b2,b3
	dvXdX,dvXdY,dvXdZ = gradient(arr, model, z)
	#dvXdX,dvXdY,dvXdZ = gradient(vec[0], model, z)
	#dvYdX,dvYdY,dvYdZ = gradient(vec[1], model, z)
	#if len(vec)==3:
	#	dvZdX,dvZdY,dvZdZ = gradient(vec[2], model, z)
	#else:
	#	dvZdX = 0.
	#	dvZdY = 0.
	#	dvZdZ = 0.
	advX = v1*dvXdX+v2*dvXdY+v3*dvXdZ
	#advY = v1*dvYdX+v2*dvYdY+v3*dvYdZ
	#advZ = v1*dvZdX+v2*dvZdY+v3*dvZdZ
	return advX
	#return (advX, advY, advZ)

def advection2(vec, model=None, dq=None, tau=None, z=None):
	'''
	Computes advection of vec.
	'''
	if not model: model=_model
	if dq == None: dq=bool(_parfile)
	if dq:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
	else: tau1=int(np.size(model.z.xc3)/2.)
	dim = len(vec)
	if dim == 3 and tau != None: to2d=True
	else: to2d=False
	if tau == None: tau=1.
	if to2d:
		v1,v2,v3 = model.z.v1[:,:,tau1],model.z.v2[:,:,tau1],model.z.v3[:,:,tau1]
	else:
		v1,v2,v3 = model.z.v1,model.z.v2,model.z.v3
	dvXdX,dvXdY,dvXdZ = gradient(vec[0], model, dq, tau, z)
	dvYdX,dvYdY,dvYdZ = gradient(vec[1], model, dq, tau, z)
	if dim == 3:
		dvZdX,dvZdY,dvZdZ = gradient(vec[2], model, dq, tau, z)
	else:
		dvZdX = 0.
		dvZdY = 0.
		dvZdZ = 0.
	advX = v1*dvXdX+v2*dvXdY+v3*dvXdZ
	advY = v1*dvYdX+v2*dvYdY+v3*dvYdZ
	advZ = v1*dvZdX+v2*dvZdY+v3*dvZdZ
	return (advX, advY, advZ)

def variation_sphere(p, model=None):
	'''
	Returns in spherical the relative position of the nMBP mesh with
	respect to the column p.

	Not used.
	'''
	if not model: model=_model
	pX, pY, pZ = model.z.xc1[p[0],0,0], model.z.xc2[0,p[1],0], model.z.xc3[0,0,p[2]]
	X, Y, Z = meshgrid(model)
	X, Y, Z = X-pX, Y-pY, Z-pZ
	r = np.sqrt(X**2+Y**2+Z**2)
	theta = np.empty_like(r)
	non_zero_r = np.where(r!=0.)
	zero_r = np.where(r==0.)
	theta[non_zero_r] = np.arccos(Z[non_zero_r]/r[non_zero_r])
	theta[zero_r] = 0.
	non_zero_x = np.where(X!=0.)
	zero_x = np.where(X==0.)
	phi = np.empty_like(r)
	phi[non_zero_x] = np.mod(.5*np.pi*(np.sign(Y[non_zero_x])-1.)+np.arctan(Y[non_zero_x]/X[non_zero_x]),2.*np.pi)
	phi[zero_x] = np.mod(.5*np.sign(Y[zero_x])*np.pi,2.*np.pi)
	return r, theta, phi

def a_centripetal(s, model=None, z=None):
	'''
	Computes centripetal acceleration from rotation around column s
	(s is given by the "snake" module).
	'''
	if not model: model=_model
	if z != None:
		pX, pY = model.z.xc1[s[0],0,0], model.z.xc2[0,s[1],0]
		X, Y = meshgrid2D()
		dX, dY = X[:,:]-pX, Y[:,:]-pY
		r2 = dX**2+dY**2
		r2[s[0],s[1]] = 1.e-10*np.min(r2[np.where(r2!=0.)])
		v2 = (-dY*model.z.v1[:,:,z]+dX*model.z.v2[:,:,z])**2/r2
	else:
		X, Y, Z = meshgrid()
		if type(s[0]) != list:
			pX, pY = model.z.xc1[s[0],0,0], model.z.xc2[0,s[1],0]
			dX, dY = X[:,:,:]-pX, Y[:,:,:]-pY
			r2 = dX**2+dY**2
			r2[s[0],s[1],:] = 1.e-10*np.min(r2[np.where(r2!=0.)])
			v2 = (-dY*model.z.v1+dX*model.z.v2)**2/r2
		else:
			dX = np.empty_like(X)
			dY = np.empty_like(Y)
			r2 = np.ones_like(X)
			v2 = np.zeros_like(Y)
			for p0 in s:
				pX, pY = model.z.xc1[p0[0],0,0], model.z.xc2[0,p0[1],0]
				dX[:,:,p0[2]], dY[:,:,p0[2]] = X[:,:,p0[2]]-pX, Y[:,:,p0[2]]-pY
				r2[:,:,p0[2]] = dX[:,:,p0[2]]**2+dY[:,:,p0[2]]**2
				r2[p0[0],p0[1],p0[2]] = 1.e-10*np.min(r2[:,:,p0[2]][np.where(r2[:,:,p0[2]]!=0.)])
				v2[:,:,p0[2]] = (-dY[:,:,p0[2]]*model.z.v1[:,:,p0[2]]+dX[:,:,p0[2]]*model.z.v2[:,:,p0[2]])**2/r2[:,:,p0[2]]
	return -v2/r2*dX, -v2/r2*dY

def a_centripetal2(model=None, s=None, tau=None, r=3, centre='rho', z=None):
	'''
	Computes centripetal acceleration from rotation around column s
	(s is given by the "snake" module).

	TODO: from 3D to 2D
	'''
	if not model: model=_model
	if z != None: to2d=True
	else: to2d=False
	if tau == None: tau=1.
	if s == None:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
		s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	if z == 'tau': z=tau1
	X, Y, Z = meshgrid(model)
	v2 = model.z.v1**2+model.z.v2**2+model.z.v3**2
	aX = np.zeros_like(v2)
	aY = np.zeros_like(v2)
	pint = np.array(s, dtype=int)
	for p0 in pint:
		if to2d and p0[2] != z: continue
		if centre == 'rho':
			pX, pY, pZ = model.z.xc1[p0[0],0,0], model.z.xc2[0,p0[1],0], model.z.xc3[0,0,p0[2]]
		else:
			pX, pY, pZ = model.z.xc1[centre[0],0,0], model.z.xc2[0,centre[1],0], model.z.xc3[0,0,p0[2]]
		dX, dY = X[:,:,p0[2]]-pX, Y[:,:,p0[2]]-pY
		r2 = dX**2+dY**2
		if centre == 'rho': r2[p0[0],p0[1]] = 1.e-10
		else: r2[centre[0],centre[1]] = 1.e-10
		aX[:,:,p0[2]] = -v2[:,:,p0[2]]/r2*dX
		aY[:,:,p0[2]] = -v2[:,:,p0[2]]/r2*dY
	if to2d: return aX[:,:,z], aY[:,:,z]
	else: return aX, aY

def a_gradient(s, model=None, z=None):
	'''
	Computes the actual acceleration due to gradient pressure
	'''
	if not model: model=_model

	grad = gradient(model.dq.P, model, z)
	grad = grad[0]**2+grad[1]**2


	if z != None:
		pX, pY = model.z.xc1[s[0],0,0], model.z.xc2[0,s[1],0]
		X, Y = meshgrid2D()
		dX, dY = X[:,:]-pX, Y[:,:]-pY
		r2 = dX**2+dY**2
		r2[s[0],s[1]] = 1.e-10*np.min(r2[np.where(r2!=0.)])
	else:
		X, Y, Z = meshgrid()
		if type(s[0]) != list:
			pX, pY = model.z.xc1[s[0],0,0], model.z.xc2[0,s[1],0]
			dX, dY = X[:,:,:]-pX, Y[:,:,:]-pY
			r2 = dX**2+dY**2
			r2[s[0],s[1],:] = 1.e-10*np.min(r2[np.where(r2!=0.)])
		else:
			dX = np.empty_like(X)
			dY = np.empty_like(Y)
			r2 = np.ones_like(X)
			for p0 in s:
				pX, pY = model.z.xc1[p0[0],0,0], model.z.xc2[0,p0[1],0]
				dX[:,:,p0[2]], dY[:,:,p0[2]] = X[:,:,p0[2]]-pX, Y[:,:,p0[2]]-pY
				r2[:,:,p0[2]] = dX[:,:,p0[2]]**2+dY[:,:,p0[2]]**2
				r2[p0[0],p0[1],p0[2]] = 1.e-10*np.min(r2[:,:,p0[2]][np.where(r2[:,:,p0[2]]!=0.)])
	r = np.sqrt(r2)
	dX, dY = dX/r, dY/r



	#grad = gradient(model.dq.P, model, z)
	if z != None: rho = model.z.rho[:,:,z]
	else: rho = model.z.rho
	return -grad*dX/rho, -grad*dY/rho
	#return -grad[0]*dX/rho, -grad[1]*dY/rho, -grad[2]/rho

def a_breaking_ratio(model=None, s=None, dq=None, odb=None, tau=1., z=None, r=3, full_advection=False, only_boundary=True, centre='rho', advect=False):
	'''
	Computes the "breaking ratio" of gradient p to centripetal
	acceleration: (a_g-a_c)/a_c.

	For an ideal cylindrically symetric field rotating around a
	given centre, this ratio is zero.

	If full_advection is True, replaces a_c by advection term.
	'''
	if not model: model=_model
	if dq == None: dq=bool(_parfile)
	if dq:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
		tau1min=int(round(np.min(tau1_level)))
	else: tau1=int(np.size(model.z.xc3)/2.)
	if s == None:
		s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	if z == None or z == 'tau': pass
	if z == 'mintau': tau1 = tau1min
	if z == 'halftau': tau1 = int(round(.5*(tau1min+tau1)))
	tau1_sindex = int(list(s[:,2]).index(tau1))
	xmean = int(s[tau1_sindex,0])
	ymean = int(s[tau1_sindex,1])
	if not odb: odb=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean))
	o, d, b = odb
	a_g = a_gradient((xmean,ymean), model,tau1)
	if not advect:
		a_c = a_centripetal((xmean,ymean),model,tau1)
		#a_c = (a_c[0][:,:,tau1], a_c[1][:,:,tau1])
		a_c = np.sqrt(a_c[0]**2+a_c[1]**2)
		#a_g = (a_g[0][:,:,tau1], a_g[1][:,:,tau1])
	else:
		a_c1 = advection(model.z.v1,model,tau1)
		a_c2 = advection(model.z.v1,model,tau1)
		a_c3 = advection(model.z.v1,model,tau1)
		a_c = np.sqrt(a_c1**2+a_c2**2+a_c3**2)
		#a_g = (a_g[0][:,:,tau1], a_g[1][:,:,tau1], a_g[2][:,:,tau1])
	a_g = np.sqrt(a_g[0]**2+a_g[1]**2)
	if only_boundary:
		return np.sqrt(np.mean(((a_g[d]-a_c[d])/a_c[d])**2)),np.max(np.abs((a_g[d]-a_c[d])/a_c[d]))
	else:
		o_star = np.zeros_like(a_c,dtype=bool)
		o_star[o] = True
		o_star[xmean,ymean] = False
		o_star=np.where(o_star)
		return np.sqrt(np.mean(((a_g[o_star]-a_c[o_star])/np.mean(a_c[o_star]))**2)),np.max(np.abs((a_g[o_star]-a_c[o_star])/np.mean(a_c[o_star])))

def wilson_depression(model=None,s=None,odb=None,tau=1.,r=3):
	'''
	Computes the Wilson depression.
	'''
	if not model: model=_model
	tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
	tau1=int(round(np.mean(tau1_level)))
	if s is None: s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	tau1_sindex = int(list(s[:,2]).index(tau1))
	xmean = int(s[tau1_sindex,0])
	ymean = int(s[tau1_sindex,1])
	if not odb: o,d,b=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean))
	else: o,d,b=odb
	min_tau1=int(np.min(tau1_level[o]))
	extended_b = np.ones_like(model.z.rho[:,:,0],dtype=bool)
	extended_b[o] = False
	extended_b = np.where(extended_b)
	mean_tau1=int(round(np.mean(tau1_level[extended_b])))
	return model.z.xc3[0,0,mean_tau1]-model.z.xc3[0,0,min_tau1]

def contrast(arr,model=None,s=None,odb=None,tau=1.,r=3):
	'''
	Computes local and global contrast of arr; (arr_i-<arr>)/<arr>.
	<arr> is taken in the local neighbourhood and on the whole z=cst
        slice, at <tau>=tau (=1 by default).
	'''
	if not model: model=_model
	tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
	tau1=int(round(np.mean(tau1_level)))
	if s is None: s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	tau1_sindex = int(list(s[:,2]).index(tau1))
	xmean = int(s[tau1_sindex,0])
	ymean = int(s[tau1_sindex,1])
	if not odb: o,d,b=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean))
	else: o,d,b=odb
	arr_mean=np.mean(arr[:,:,tau1])
	arr_bmean=np.mean(arr[:,:,tau1][b])
	contrast_global = (arr[xmean,ymean,tau1]-arr_mean)/arr_mean
	contrast_local = (arr[xmean,ymean,tau1]-arr_bmean)/arr_bmean
	return contrast_global, contrast_local

def contrastf(arr,odb,f,fOnWholeArr=True):
	'''
	Computes local and global contrast of arr; (f(arr)-<arr>)/<arr>.
	<arr> is taken in the local neighbourhood.
	'''
	o,d,b     = odb
	if fOnWholeArr:
		arr_max   = f(arr)
	else:
		arr_max   = f(arr[o])
	arr_mean  = np.mean(arr)
	arr_bmean = np.mean(arr[b])
	contrast_global = (arr_max-arr_mean)/arr_mean
	contrast_local  = (arr_max-arr_bmean)/arr_bmean
	return contrast_global, contrast_local

def refine(arr, n):
	'''
	This function interpolates an array on a finer grid
	'''
	assert(type(n) is int)
	s = np.shape(arr)
	dim = len(s)
	assert(dim==1 or dim==2)
	if dim==2:
		px, py = np.where(np.ones_like(arr))
		data = arr[px, py]
		w = (n*px, n*py)
		grid_x, grid_y = np.mgrid[:n*(s[0]-1)+1,:n*(s[1]-1)+1]
		return griddata(w, data, (grid_x, grid_y), method='linear')
	if dim==1:
		x = np.arange(s[0])
		f = interp1d(x,arr)
		return f(np.arange(n*(s[0]-1)+1)/float(n))


def plotv_slice(arr,model=None,s=None,dq=None,tau=1.,r=3,boxtext=None,show=True,rf=10,tight=True):
	'''
	Plots vertical slices of arr at y=cst and x=cst respectively.
	Depending on optional parameters, also plots tau=tau (=1 by
	default) isosurface. Marks in red the "eye" of the nMBP.
	'''
	if not model: model=_model
	if dq == None: dq=bool(_parfile)
	if dq:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
	else: tau1=int(np.size(model.z.xc3)/2.)
	if s is None: s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	sx, sy, sz = np.shape(arr)
	x, y, z = np.array(s, dtype=int).T
	xc, yc, zc = model.z.xc1[:,0,0]/1.e5, model.z.xc2[0,:,0]/1.e5, model.z.xc3[0,0,:]/1.e5
	if dq:
		tau1_sindex = int(list(s[:,2]).index(tau1))
		xmean = int(s[tau1_sindex,0])
		ymean = int(s[tau1_sindex,1])
	else:
		xmean = int(round(np.mean(x)))
		ymean = int(round(np.mean(y)))
		tau1 = 0
	sharey = False
	if tight: z0 = z
	else:
		if dq:
			xp,yp,z1p,z2p = xc, yc, zc[tau1_level[:,ymean]]-zc[tau1],zc[tau1_level[xmean,:]]-zc[tau1]
			z0min = z[0]
			z0max = z[-1]+1
			z0min = min(z0min,np.min(tau1_level[:,ymean]))
			z0max = max(z0max,1+np.max(tau1_level[:,ymean]))
			z0min = min(z0min,np.min(tau1_level[xmean,:]))
			z0max = max(z0max,1+np.max(tau1_level[ymean,:]))
			delta = z0max - z0min
			z0min = max(0, int(z0min-.1*delta))
			z0max = min(sz, int(z0max+.1*delta))
			z0 = np.arange(z0min, z0max)
			sharey=True
		else:
			z0 = range(sz)
	sliceY=arr[np.ix_(range(sx),[ymean],z0)][:,0,:].T
	sliceX=arr[np.ix_([xmean],range(sy),z0)][0,:,:].T
	vY=model.z.v2[np.ix_(range(sx),[ymean],z0)][:,0,:].T
	vX=model.z.v1[np.ix_([xmean],range(sy),z0)][0,:,:].T
	if rf > 1:	# Refine factor
		vY = refine(vY,8)
		vX = refine(vX,8)
	vY=(vY>=0.)
	vX=(vX>=0.)
	extY=np.array([xc[0],xc[sx-1],zc[z0[0]]-zc[tau1],zc[z0[-1]]-zc[tau1]])
	extX=np.array([yc[0],yc[sy-1],zc[z0[0]]-zc[tau1],zc[z0[-1]]-zc[tau1]])
	f, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,5), sharey=sharey)
	ax1.set_aspect('equal')
	ax2.set_aspect('equal')
	plt.tight_layout(pad=4)
	if _modelfile:
		title = _modelfile.split('.')
		title = '.'.join(title[:max(1,len(title)-1)])
		plt.suptitle(title,fontsize=16)
	vmin = min(np.min(sliceY), np.min(sliceX))
	vmax = max(np.max(sliceY), np.max(sliceX))
	im1 = ax1.imshow(sliceY, origin='bottom',alpha=1.,cmap=None, extent=extY, vmin=vmin, vmax=vmax)
	im2 = ax2.imshow(sliceX, origin='bottom',alpha=1.,cmap=None, extent=extX, vmin=vmin, vmax=vmax)
	f.subplots_adjust(right=0.8)
	cbar_ax = f.add_axes([0.83, 0.15, 0.02, 0.69])
	cbar = f.colorbar(im2, cax=cbar_ax)
	cbar.set_label('Density [g$\,$cm$^{-3}$]', rotation=270, labelpad=20)
	ax1.imshow(-vY, origin='bottom',alpha=.4,cmap=cm.gray, extent=extY)
	ax2.imshow(vX, origin='bottom',alpha=.4,cmap=cm.gray, extent=extX)
	xp,yp,zp = xc[x],yc[y],zc[z]-zc[tau1]
	if rf > 1:
		xp,yp,zp = refine(xp,rf), refine(yp,rf), refine(zp,rf)
		xp = savitzky_golay(xp, rf**2+1, 3)
		yp = savitzky_golay(yp, rf**2+1, 3)
		zp = savitzky_golay(zp, rf**2+1, 3)
	ax1.plot(xp,zp,'r')
	ax2.plot(yp,zp,'r')
	#ax1.plot(xc[x],zc[z]-zc[tau1],'r')
	#ax2.plot(yc[y],zc[z]-zc[tau1],'r')
	if dq:
		xp,yp,z1p,z2p = xc, yc, zc[tau1_level[:,ymean]]-zc[tau1],zc[tau1_level[xmean,:]]-zc[tau1]
		if rf > 1:
			xp,yp,z1p,z2p = refine(xp,rf), refine(yp,rf), refine(z1p,rf), refine(z2p,rf)
			xp = savitzky_golay(xp, rf*(rf-1)+1, 3)
			yp = savitzky_golay(yp, rf*(rf-1)+1, 3)
			z1p = savitzky_golay(z1p, rf*(rf-1)+1, 3)
			z2p = savitzky_golay(z2p, rf*(rf-1)+1, 3)
		ax1.plot(xp,z1p,'b')
		ax2.plot(yp,z2p,'b')
	#	ax1.plot(xc,zc[tau1_level[:,ymean]]-zc[tau1],'b')
	#	ax2.plot(yc,zc[tau1_level[xmean,:]]-zc[tau1],'b')
	ax1.set_xlim((extY[0],extY[1]))
	ax1.set_ylim((extY[2],extY[3]))
	ax2.set_xlim((extX[0],extX[1]))
	ax2.set_ylim((extX[2],extX[3]))
	ax1.set_xlabel("Spatial horizontal position (X axis) [km]")
	ax2.set_xlabel("Spatial horizontal position (Y axis) [km]")
	ax1.set_ylabel("Height [km]")
	if boxtext:
		props = dict(boxstyle='round, pad=.7, rounding_size=.2', facecolor='white', edgecolor='black', alpha=.8)
		ax1.text(0.05, 0.9, boxtext, size='smaller', ha='left', va='top', transform=ax1.transAxes, bbox=props)
	if show: plt.show()
	return sliceY, sliceX

def report(ireport=None, model=None, s=None, dq=None, odb=None, tau=1., z=None, r=3, savefig=False, show=True):
	if not model: model=_model
	if dq == None: dq=bool(_parfile)
	if dq:
		tau1_level=np.array(pybold.level(model.dq.tau,tau),dtype=int)
		tau1=int(round(np.mean(tau1_level)))
		tau1min=int(round(np.min(tau1_level)))
	else: tau1=int(np.size(model.z.xc3)/2.)
	if s == None: s=snake.snake_from_box(model.z.rho,radius=r,start=tau1)
	if z == None or z == 'tau': pass
	if z == 'mintau': tau1 = tau1min
	if z == 'halftau': tau1 = int(round(.5*(tau1min+tau1)))
	tau1_sindex = int(list(s[:,2]).index(tau1))
	xmean = int(s[tau1_sindex,0])
	ymean = int(s[tau1_sindex,1])
	print(str(xmean)+', '+str(ymean))
	print(str(np.unravel_index(np.argmin(model.z.v1[46:56,39:50,tau1]**2+model.z.v2[46:56,39:50,tau1]**2),(10,11))))
	if not odb: odb=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean),threshold=.5)
	rot_index = np.argmin(model.z.v1[:,:,tau1][odb[0]]**2)
	xrot,yrot = odb[0][0][rot_index], odb[0][1][rot_index]
	print(str(xrot)+', '+str(yrot))
	#
	#if not odb: odb=snake.select_disk(tau1_level,(xmean,ymean),mask_flag=False,threshold=.8)
	plt.imshow(model.z.rho[:,:,tau1],origin='bottom')
	oring=np.zeros_like(model.z.rho[:,:,tau1],dtype=bool)
	oring[odb[0]]=True
	plt.imshow(oring,origin='bottom',alpha=.7,interpolation='none')
	plt.quiver(model.z.v2[:,:,tau1],model.z.v1[:,:,tau1])
	#odb2=snake.select_disk(model.z.rho[:,:,tau1],(xmean,ymean),threshold=.8)
	#plt.figure()
	#plt.imshow(model.z.rho[:,:,tau1],origin='bottom')
	#oring=np.zeros_like(model.z.rho[:,:,tau1],dtype=bool)
	#oring[odb2[0]]=True
	#plt.imshow(oring,origin='bottom',alpha=.7,interpolation='none')
	#plt.quiver(model.z.v2[:,:,tau1],model.z.v1[:,:,tau1])
	#plt.figure()
	#plt.imshow(model.z.v2[:,:,tau1]**2+model.z.v1[:,:,tau1]**2,origin='bottom',interpolation='none',cmap=cm.gray)
	#
	fmt='{0: <20}'
	if _modelfile: print(fmt.format('File:')+_modelfile)
	else: print(fmt.format('File:')+'Unknown')
	wd=wilson_depression(model,s,odb,tau,r)
	print(fmt.format('Wilson depression:')+str(wd))
	#c_rho=contrast(model.z.rho,model,s,odb,tau,r)
	c_rho=contrastf(model.z.rho[:,:,tau1],odb,np.min)
	print(fmt.format('Density contrast:')+str(c_rho))
	#c_p=contrast(model.dq.P,model,s,odb,tau,r)
	c_p=contrastf(model.dq.P[:,:,tau1],odb,np.min)
	print(fmt.format('Pressure contrast:')+str(c_p))
	a_ratio=a_breaking_ratio(model,s,dq,odb,tau,z,r,centre='rho')
	print(fmt.format('Acceleration ratio:')+str(a_ratio))
	adv_ratio=a_breaking_ratio(model,s,dq,odb,tau,z,r,centre='rho',advect=True)
	print(fmt.format('Advection ratio::')+str(adv_ratio))
	if _modelfile:
		title = _modelfile.split('.')
		title = '.'.join(title[:max(1,len(title)-1)])
	if ireport and _modelfile:
		id_name = title.split('nMBP')[-1]
		with open(ireport) as f:
			lines = f.readlines()
		names = [l.split('\t')[0] for l in lines]
		lines = [[f.strip() for f in l.split('\t')] for l in lines[1:]]
		idx = int(names.index(id_name))+1
		diam = float(lines[idx][1])
		c_I = (float(lines[idx][2]),float(lines[idx][3]))
		time = float(lines[idx][4])
		print(fmt.format('Diametre:')+str(diam))
		print(fmt.format('Intensity contrast:')+str(c_I))
		print(fmt.format('Lifetime:')+str(time))
	boxtext = 'Wilson depression: '+str(int(round(wd/1.e5)))+' [km]\n'
	boxtext+= 'Density contrast: '+str(int(round(100*c_rho[0])))+'%, '
	boxtext+= str(int(round(100*c_rho[1])))+'%\n'
	boxtext+= 'Pressure contrast: '+str(int(round(100*c_p[0])))+'%, '
	boxtext+= str(int(round(100*c_p[1])))+'%'
	#boxtext+= 'Acceleration ratio: '+str(int(round(100*a_ratio[0])))+'%, '
	#boxtext+= str(int(round(100*a_ratio[1])))+'%'
	if ireport and _modelfile:
		delta = model.z.xc1[1,:,:]-model.z.xc1[0,:,:]
		boxtext+= '\nDiametre: '+str(int(round(diam*delta/1.e5)))+' [km]'
		boxtext+= '\nIntensity contrast: '+str(int(round(100*c_I[0])))+'%, '
		boxtext+= str(int(round(100*c_I[1])))+'%'
		boxtext+= '\nLifetime: '+str(int(round(time)))+' [s]'
	if savefig:
		sY,sX=plotv_slice(model.z.rho,model,s,dq,tau,r,boxtext,show=False)
		if savefig == True: plt.savefig(title+'.png')
		else: plt.savefig(savefig)
		if show:
			plt.show()
		plt.close()
	else:
		sY,sX=plotv_slice(model.z.rho,model,s,dq,tau,r,boxtext,show=show)
		plt.close()
