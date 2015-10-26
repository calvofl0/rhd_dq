#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python interface to the rhd_dq backend Fortran module.

Getting started
===============

To load the third snapshot of 'rhd.full' you may proceed as follows (we
assume parfile is 'rhd.par' and contains the correct path to EOS and OPA
files):

>>> model = uio_struct()
>>> model.load('rhd.full', 'rhd.par', imodel=3)

Now you can print a map with the contents:

>>> model
dataset_id          <type 'str'>
dq                  <class 'pybold.uio_struct'>
dtime               <type 'float'>
head                <class 'pybold.uio_struct'>
modelitime          <type 'int'>
modeltime           <type 'float'>
time_out_full_last  <type 'float'>
time_out_mean_last  <type 'float'>
z                   <class 'pybold.uio_struct'>

>>> model.z
box_id              <type 'str'>
dimension           <type 'numpy.ndarray'>        (2, 3)
ei                  <type 'numpy.ndarray'>        (140, 140, 150)
itime               <type 'int'>
rho                 <type 'numpy.ndarray'>        (140, 140, 150)
time                <type 'float'>
v1                  <type 'numpy.ndarray'>        (140, 140, 150)
v2                  <type 'numpy.ndarray'>        (140, 140, 150)
v3                  <type 'numpy.ndarray'>        (140, 140, 150)
xb1                 <type 'numpy.ndarray'>        (141, 1, 1)
xb2                 <type 'numpy.ndarray'>        (1, 141, 1)
xb3                 <type 'numpy.ndarray'>        (1, 1, 151)
xc1                 <type 'numpy.ndarray'>        (140, 1, 1)
xc2                 <type 'numpy.ndarray'>        (1, 140, 1)
xc3                 <type 'numpy.ndarray'>        (1, 1, 150)

To display available derived quantities, you can type:

>>> model.dq
B2                  <type 'NoneType'>
Babs                <type 'NoneType'>
Bc1                 <type 'NoneType'>
Bc2                 <type 'NoneType'>
Bc3                 <type 'NoneType'>
Bh                  <type 'NoneType'>
P                   <type 'NoneType'>
T                   <type 'NoneType'>
beta                <type 'NoneType'>
...

What is P?

>>> model.dq['P'].desc
'Pressure'

We can also change the display mode:

>>> model.names
>>> model.dq.P
'Pressure'

Let's come back to the default 'values' mode:

>>> model.values

Let's compute pressure:

>>> model.dq.P
 P         [ok]
array([[[  4.81998250e+06,   4.67989200e+06,   4.53605400e+06, ...,
           2.92639351e+01,   2.44200096e+01,   2.03071632e+01],
        [  4.81881050e+06,   4.67870050e+06,   4.53405600e+06, ...,
           2.84149323e+01,   2.35653324e+01,   1.93235435e+01],
        [  4.81854850e+06,   4.67866950e+06,   4.53386150e+06, ...,
           2.71816387e+01,   2.23898506e+01,   1.80947742e+01],
        ...,
        [  4.82060800e+06,   4.67813850e+06,   4.53572550e+06, ...,
           3.15569439e+01,   2.66288719e+01,   2.26559296e+01],
        [  4.82066750e+06,   4.67926750e+06,   4.53648250e+06, ...,
           3.06190186e+01,   2.57513676e+01,   2.17916222e+01],
        [  4.82049750e+06,   4.67987900e+06,   4.53652150e+06, ...,
           2.99480400e+01,   2.51007214e+01,   2.10807152e+01]],
	...
        [  4.81750100e+06,   4.67527500e+06,   4.53174750e+06, ...,
           3.02413216e+01,   2.53970299e+01,   2.13400726e+01]]], dtype=float32)

Pressure is now in memory and will not be re-computed if requested again.

We can now print the units for pressure:

>>> model.dq['P'].unit
'erg/g*g/cm^3'

This show that pressure is energy per mass times density, which is energy
per volume. If we were in the units mode:

>>> model.units

we could directly print units with

>>> model.dq.P
'erg/g*g/cm^3'

Let's open now a mean file. We can either create a new uio_struct instance,
or we can empty the model variable:

>>> model.record('delete')

Now we load 'rhd.mean':

>>> model.load('rhd.mean')

Note that 'rhd.par' has not been provided this time because we do not
wish to compute derived quantities (indeed, derived quantities can only
be computed for FULL files). The rhd_dq backend will therefore not be
used.

Important notes concerning rhd_dq backend
=========================================

- The rhd_dq backend requires an important stack size to run, otherwise
  a segmentation fault is to be expected. In linux one can use the
  command 'ulimit -s unlimited' before starting python to fix this.

- Only one instance of the backend is running at a time. This means that
  one can run multiple instances of uio_struct as long a only one of them
  uses rhd_dq.

- How do we now when an uio_struct instance is using rhd_dq? Providing a
  parfile to the load or the open function will automatically initialize
  the rhd_dq backend, otherwise it will not be used.

- When rhd_dq backend is used, magnetic fields are expressed in Gauss.
  /!\ However, without rhd_dq backend there is no post-processing at all
  after reading models, and magnetic field units are therefore
  G/sqrt(4*pi)
"""

import numpy as np
import sys
from _pybold import *
try:
	from pypmd import *
	with_pypmd=True
except ImportError:
	with_pypmd=False

def namespace():
	'''
	Create an empty module to be used as generic namespace.

	>>> ns = namespace()
	>>> dir(ns)
	[]

	>>> ns.member = 'value'
	>>> ns.member
	'value'
	'''

	import imp
	ns = imp.new_module('')
	for attr in dir(ns):
		delattr(ns, attr)
	return ns

# Mathematical and Physical constants in cgs units
pi	= np.pi
sqrt4pi	= np.sqrt(4*np.pi)
hbar	= 1.054571726e-27
c	= 2.99792458e10
kB	= 1.3806488e-16
sigma	= pi**2*kB**4/(60*hbar**3*c**2)
a	= 4*sigma/c
NA	= 6.02214129e23
R	= kB*NA
alpha	= 7.2973525698e-3
e	= np.sqrt(hbar*alpha*c)
G	= 6.67384e-8
amu	= 1.660538921e-24

# SI conversion
si	= namespace()
si.cm	= 1.e-2
si.g	= 1.e-3
si.s	= 1.
si.dyne	= 1.e-5
si.erg	= 1.e-7
si.barye= 1.e-1
si.m0	= 4.*pi*1.e-7
si.e0	= 1./(si.m0*(c*(si.cm/si.s))**2)
si.e	= e*np.sqrt(si.erg*si.cm*4.*pi*si.e0)

# Natural units conversion
nu	= namespace()
nu.erg	= si.erg/si.e
nu.cm	= 1/(hbar*c*nu.erg)
nu.s	= c*nu.cm
nu.g	= c**2*nu.erg

# Solar constants (in cgs units)
sun	= namespace()
sun.M	= 1.98855e33
sun.R 	= 6.96342e10
sun.g	= 27542.29
sun.T	= 5778
sun.Tc	= 1.57e7
sun.L	= 3.846e31
sun.I	= sun.L/(2*pi*sun.R)**3
sun.rhoc= 1.622e2
sun.rhom= 1.408
sun.D	= 1.496e13
sun.mu	= 0.594

# Level function
# Given an array arr, some axis ax and a level l, returns an array
# privated of axis ax, whose elements are the coordinate on ax at which
# arr is equal to l
# Every column parametrized along axis ax must be strictly monotonic

def level1D(col, l, t=None, rounding=False):
	"""
	Given a monotonic vector 'col', returns the index corresponding
	to value 'l' when t is an integer type.

	If t is a floating point type, a linear interpolation is computed
	to return a "floating point index" of level l.

	When type is integer, if rounding is True, the chosen index is the
	closest to level l (typically level l is between two indices).
	"""
	if not t: t=np.array(col, copy=False).dtype.type
	ext=(col<=l)
	ext=np.logical_or(ext[:-1], ext[1:])
	ext=np.where(ext)
	if len(ext)!=1: return np.nan
	ext=ext[0][0]
	if rounding:
		return t(ext)+t(round((t(l)-col[ext])/(col[ext+1]-col[ext])))
	return t(ext)+t((t(l)-col[ext])/(col[ext+1]-col[ext]))

def level(arr, l, ax=2, t=None, rounding=False): 
	"""
	As level 1D but accepts general arrays and level is taken is some
	specified axis.
	"""
	return np.apply_along_axis(level1D, ax, arr, l, t, rounding)

def varAtMeanLevel(var, arr, l, ax=2):
	"""
	Returns the values of var in axis ax corresponding to the position in
	which arr is in average of level l.
	"""
	isosurf=level(arr,l,ax)
	average=int(round(np.mean(isosurf)))
	average=average*np.ones_like(isosurf, dtype=int)
	nx,ny=np.shape(isosurf)
	indices=list(np.meshgrid(range(nx),range(ny)))+[average]
	lookup=[1,0]
	lookup.insert(ax,2)
	indices=tuple([indices[i] for i in lookup])
	return var[indices]

def varAtLevel(var, arr, l, ax=2, interpolate=True):
	"""
	Returns the values of var in axis ax corresponding to the position in
	which arr is of level l, using the level function to compute
	isosurface.

	Values of var are linearly interpolated unless interpolate is set to
	False.
	"""
	if not interpolate:
		#isosurf=level(arr,l,ax,t=int,rounding=True)
		isosurf=level(arr,l,ax).astype(int)
	else:
		fisosurf=level(arr,l,ax)
		frac=fisosurf%1.
		isosurf=(fisosurf-frac).astype(int)
	nx,ny=np.shape(isosurf)
	indices=list(np.meshgrid(range(nx),range(ny)))+[isosurf]
	if interpolate:
		indicesR=list(np.meshgrid(range(nx),range(ny)))+[isosurf+1]
	lookup=[1,0]
	lookup.insert(ax,2)
	indices=tuple([indices[i] for i in lookup])
	if not interpolate:
		return var[indices]
	else:
		indicesR=tuple([indicesR[i] for i in lookup])
		return (1.-frac)*var[indices]+frac*var[indicesR]

def local_minima(arr, threshold=.2):
	"""
	Returns indices of local minimal of a 2D array
	"""
	import scipy.ndimage as ndimage
	import scipy.ndimage.filters as filters
	data_max=(filters.maximum_filter(arr, size=3)==arr)
	# To be continued...

def extend_periodically(vector, period, direction):
	s = np.shape(vector)
	flatvect = vector.flatten()
	l = len(flatvect)
	out = np.array(flatvect, copy=True, order='F')
	if direction>0:
		for i in range(l-1):
			if flatvect[i+1]<flatvect[i]:
				out[i+1:]+=period
				break
	elif direction<0:
		for i in range(l-1,0,-1):
			if flatvect[i-1]>flatvect[i]:
				out[:i]-=period
				break
	return np.asfortranarray(out.reshape(s))

class uio_struct_item(object):
	__slots__ = ['name', 'desc', 'unit', 'value', 'self', '__link', 'record', '_parent']
	@property
	def self(self):
		return self
	def __init__(self, record=None):
		object.__setattr__(self, 'name', None)
		object.__setattr__(self, 'desc', None)
		object.__setattr__(self, 'unit', None)
		object.__setattr__(self, 'value', None)
		object.__setattr__(self, '_uio_struct_item__link', False)
		object.__setattr__(self, 'record', record)
	def get_dq(self, dq) :
		"""
	Get derived quantity from Fortran rhd_dq backend. Unless you are
	debugging rhd_dq backend module, you are not intended to use
	this function.
		"""
		link = self._uio_struct_item__link
		self._uio_struct_item__link = False
		if self.value == 'linked':
			self._uio_struct_item__link = link
		else : get_dq(self.record, dq)
	def __getattribute__(self, name):
		if object.__getattribute__(self, '_uio_struct_item__link') and name=='value':
			if self.name == 'rho' : return box.dq_rho
			elif self.name == 'ei' : return box.dq_ei
			elif self.name == 'v1' : return box.dq_v1
			elif self.name == 'v2' : return box.dq_v2
			elif self.name == 'v3' : return box.dq_v3
			elif self.name == 'bb1' : return box.dq_Bb1
			elif self.name == 'bb2' : return box.dq_Bb2
			elif self.name == 'bb3' : return box.dq_Bb3
			elif self.name == 'xb1' :
				return box.dq_xb1.reshape((len(box.dq_xb1),1,1))
			elif self.name == 'xb2' :
				return box.dq_xb2.reshape((1,len(box.dq_xb2),1))
			elif self.name == 'xb3' :
				return box.dq_xb3.reshape((1,1,len(box.dq_xb3)))
			elif self.name == 'xc1' :
				return box.dq_xc1.reshape((len(box.dq_xc1),1,1))
			elif self.name == 'xc2' :
				return box.dq_xc2.reshape((1,len(box.dq_xc2),1))
			elif self.name == 'xc3' :
				return box.dq_xc3.reshape((1,1,len(box.dq_xc3)))
			elif self.name == 'dimension' :
				return np.array([[box.dq_m1,box.dq_m2,box.dq_m3],[box.dq_n1,box.dq_n2,box.dq_n3]])
			elif self.name == 'time' : return box.dq_time.item()
			elif self.name == 'itime' : return box.dq_itime.item()
			elif self.name == 'T' :
				self.get_dq('T')
				return box.dq_T
			elif self.name == 'P' :
				self.get_dq('P')
				return box.dq_P
			elif self.name == 's' :
				self.get_dq('s')
				return box.dq_s
			elif self.name == 'j1' :
				self.get_dq('j1')
				return box.dq_j1
			elif self.name == 'j2' :
				self.get_dq('j2')
				return box.dq_j2
			elif self.name == 'j3' :
				self.get_dq('j3')
				return box.dq_j3
			elif self.name == 'jabs' :
				self.get_dq('jabs')
				return box.dq_jabs
			elif self.name == 'kappa' :
				self.get_dq('kappa')
				return box.dq_kappa
			elif self.name == 'Bc1' :
				self.get_dq('Bc1')
				return box.dq_Bc1
			elif self.name == 'Bc2' :
				self.get_dq('Bc2')
				return box.dq_Bc2
			elif self.name == 'Bc3' :
				self.get_dq('Bc3')
				return box.dq_Bc3
			elif self.name == 'divB' :
				self.get_dq('divB')
				return box.dq_divB
			elif self.name == 'vabs' :
				self.get_dq('vabs')
				return box.dq_vabs
			elif self.name == 'vh' :
				self.get_dq('vh')
				return box.dq_vh
			elif self.name == 'ekin' :
				self.get_dq('ekin')
				return box.dq_ekin
			elif self.name == 'plin' :
				self.get_dq('plin')
				return box.dq_plin
			elif self.name == 'vmflux' :
				self.get_dq('vmflux')
				return box.dq_vmflux
			elif self.name == 'g1' :
				self.get_dq('g1')
				return box.dq_g1
			elif self.name == 'g3' :
				self.get_dq('g3')
				return box.dq_g3
			elif self.name == 'cs' :
				self.get_dq('cs')
				return box.dq_cs
			elif self.name == 'mach' :
				self.get_dq('mach')
				return box.dq_mach
			elif self.name == 'mmu' :
				self.get_dq('mmu')
				return box.dq_mmu
			elif self.name == 'tau' :
				self.get_dq('tau')
				return box.dq_tau
			elif self.name == 'Babs' :
				self.get_dq('Babs')
				return box.dq_Babs
			elif self.name == 'Bh' :
				self.get_dq('Bh')
				return box.dq_Bh
			elif self.name == 'B2' :
				self.get_dq('B2')
				return box.dq_B2
			elif self.name == 'emag' :
				self.get_dq('emag')
				return box.dq_emag
			elif self.name == 'cA' :
				self.get_dq('cA')
				return box.dq_cA
			elif self.name == 'beta' :
				self.get_dq('beta')
				return box.dq_beta
			elif self.name == 'csca' :
				self.get_dq('csca')
				return box.dq_csca
			else : return None
			#return object.__getattribute__(self, self.__link)
		else :
			return object.__getattribute__(self, name)
	def __setattr__(self, name, value):
		def safecheck(var, value):
			if(np.shape(var) != np.shape(value) or type(var) != type(value)) :
				raise ValueError('Variable type must be '+str(type(var))+' and shape must be '+str(np.shape(var))+'.')
		if object.__getattribute__(self, '_uio_struct_item__link') and name=='value':
			if self.name == 'rho' :
				safecheck(box.dq_rho,value)
				box.dq_rho=value
			elif self.name == 'ei' :
				safecheck(box.dq_ei,value)
				box.dq_ei=value
			elif self.name == 'v1' :
				safecheck(box.dq_v1,value)
				box.dq_v1=value
			elif self.name == 'v2' :
				safecheck(box.dq_v2,value)
				box.dq_v2=value
			elif self.name == 'v3' :
				safecheck(box.dq_v3,value)
				box.dq_v3=value
			elif self.name == 'bb1' :
				safecheck(box.dq_Bb1,value)
				box.dq_Bb1=value
			elif self.name == 'bb2' :
				safecheck(box.dq_Bb2,value)
				box.dq_Bb2=value
			elif self.name == 'bb3' :
				safecheck(box.dq_Bb3,value)
				box.dq_Bb3=value
			elif self.name == 'xb1' :
				safecheck(box.dq_xb1.reshape((len(box.dq_xb1),1,1)),value)
				box.dq_xb1=value[:,0,0]
			elif self.name == 'xb2' :
				safecheck(box.dq_xb2.reshape(1,(len(box.dq_xb2),1)),value)
				box.dq_xb2=value[0,:,0]
			elif self.name == 'xb3' :
				safecheck(box.dq_xb3.reshape(1,1,(len(box.dq_xb3))),value)
				box.dq_xb3=value[0,0,:]
			elif self.name == 'xc1' :
				safecheck(box.dq_xc1.reshape((len(box.dq_xc1),1,1)),value)
				box.dq_xc1=value[:,0,0]
			elif self.name == 'xc2' :
				safecheck(box.dq_xc2.reshape(1,(len(box.dq_xc2),1)),value)
				box.dq_xc2=value[0,:,0]
			elif self.name == 'xc3' :
				safecheck(box.dq_xc3.reshape(1,1,(len(box.dq_xc3))),value)
				box.dq_xc3=value[0,0,:]
			elif self.name == 'dimension' :
				safecheck(np.empty(2,3),value)
				box.dq_m1=value[0,0]
				box.dq_n1=value[1,0]
				box.dq_m2=value[0,1]
				box.dq_n2=value[1,1]
				box.dq_m3=value[0,2]
				box.dq_n3=value[1,2]
			elif self.name == 'time' :
				box.dq_time = value
			elif self.name == 'itime' :
				box.dq_itime = value
			elif self.name == 'T' :
				self.get_dq('T')
				safecheck(box.dq_T,value)
				box.dq_T=value
			elif self.name == 'P' :
				self.get_dq('P')
				safecheck(box.dq_P,value)
				box.dq_P=value
			elif self.name == 's' :
				self.get_dq('s')
				safecheck(box.dq_s,value)
				box.dq_s=value
			elif self.name == 'j1' :
				self.get_dq('j1')
				safecheck(box.dq_j1,value)
				box.dq_j1=value
			elif self.name == 'j2' :
				self.get_dq('j2')
				safecheck(box.dq_j2,value)
				box.dq_j2=value
			elif self.name == 'j3' :
				self.get_dq('j3')
				safecheck(box.dq_j3,value)
				box.dq_j3=value
			elif self.name == 'jabs' :
				self.get_dq('jabs')
				safecheck(box.dq_jabs,value)
				box.dq_jabs=value
			elif self.name == 'kappa' :
				self.get_dq('kappa')
				safecheck(box.dq_kappa,value)
				box.dq_kappa=value
			elif self.name == 'Bc1' :
				self.get_dq('Bc1')
				safecheck(box.dq_Bc1,value)
				box.dq_Bc1=value
			elif self.name == 'Bc2' :
				self.get_dq('Bc2')
				safecheck(box.dq_Bc2,value)
				box.dq_Bc2=value
			elif self.name == 'Bc3' :
				self.get_dq('Bc3')
				safecheck(box.dq_Bc3,value)
				box.dq_Bc3=value
			elif self.name == 'divB' :
				self.get_dq('divB')
				safecheck(box.dq_divB,value)
				box.dq_divB=value
			elif self.name == 'vabs' :
				self.get_dq('vabs')
				safecheck(box.dq_vabs,value)
				box.dq_vabs=value
			elif self.name == 'vh' :
				self.get_dq('vh')
				safecheck(box.dq_vh,value)
				box.dq_vh=value
			elif self.name == 'ekin' :
				self.get_dq('ekin')
				safecheck(box.dq_ekin,value)
				box.dq_ekin=value
			elif self.name == 'plin' :
				self.get_dq('plin')
				safecheck(box.dq_plin,value)
				box.dq_plin=value
			elif self.name == 'vmflux' :
				self.get_dq('vmflux')
				safecheck(box.dq_vmflux,value)
				box.dq_vmflux=value
			elif self.name == 'g1' :
				self.get_dq('g1')
				safecheck(box.dq_g1,value)
				box.dq_g1=value
			elif self.name == 'g3' :
				self.get_dq('g3')
				safecheck(box.dq_g3,value)
				box.dq_g3=value
			elif self.name == 'cs' :
				self.get_dq('cs')
				safecheck(box.dq_cs,value)
				box.dq_cs=value
			elif self.name == 'mach' :
				self.get_dq('mach')
				safecheck(box.dq_mach,value)
				box.dq_mach=value
			elif self.name == 'mmu' :
				self.get_dq('mmu')
				safecheck(box.dq_mmu,value)
				box.dq_mmu=value
			elif self.name == 'tau' :
				self.get_dq('tau')
				safecheck(box.dq_tau,value)
				box.dq_tau=value
			elif self.name == 'Babs' :
				self.get_dq('Babs')
				safecheck(box.dq_Babs,value)
				box.dq_Babs=value
			elif self.name == 'Bh' :
				self.get_dq('Bh')
				safecheck(box.dq_Bh,value)
				box.dq_Bh=value
			elif self.name == 'B2' :
				self.get_dq('B2')
				safecheck(box.dq_B2,value)
				box.dq_B2=value
			elif self.name == 'emag' :
				self.get_dq('emag')
				safecheck(box.dq_emag,value)
				box.dq_emag=value
			elif self.name == 'cA' :
				self.get_dq('cA')
				safecheck(box.dq_cA,value)
				box.dq_cA=value
			elif self.name == 'beta' :
				self.get_dq('beta')
				safecheck(box.dq_beta,value)
				box.dq_beta=value
			elif self.name == 'csca' :
				self.get_dq('csca')
				safecheck(box.dq_csca,value)
				box.dq_csca=value
			#object.__setattr__(self, self.__link, value)
		else :
			object.__setattr__(self, name, value)

	#def __getattribute__(self, name):
	#	print(name)
	#	if hasattr(self, name) :
	#		print('has')
	#		return 'End' #object.__getattribute__(self, name)
	#	else :
	#		return 'End2' #object.__getattribute__(self.value, name)
	#def __setattr__(self, name, value):
	#	if hasattr(self, name):
	#		object.__setattr__(self, name, value)
	#		self.__dict__[name] = value
	#	else :
	#			raise AttributeError("'"+type(self.value).__name__+"' object has no attribute '"+name+"'")
	def __repr__(self):
		return(self.value.__repr__())

class uio_table(list):
	__slots__ = ['headers', 'desc', 'units', 'table']
	def __init__(self):
		object.__setattr__(self, 'headers', [])
		object.__setattr__(self, 'desc', [])
		object.__setattr__(self, 'units', [])
		object.__setattr__(self, 'table', [])
	#def __setattr__(self, name, value):
		#if hasattr(self, name):
		#	object.__setattr__(self, name, value)
		#	self.__dict__[name] = value
		#else :
		#		raise AttributeError("'"+type(self.value).__name__+"' object has no attribute '"+name+"'")
	def __repr__(self):
		istable = (self.table.__class__ == list)
		if istable :
			for col in self.table:
				if col.__class__!=list and col.__class__!=np.ndarray:
					istable = False
		if istable:
			_repr_ = ''
			if self.headers.__class__ == list:
				if len(self.headers)>0:
					for i in range(len(self.headers)):
						_repr_ += str(self.headers[i])
						_repr_ += '\t'
					_repr_ += '\n'
			if len(self.table) > 0:
				sz = max([len(self.table[i]) for i in range(len(self.table))])
			else: sz = 0
			for i in range(sz):
				for j in range(len(self.table)):
					if i < len(self.table[j]):
						_repr_ += str(self.table[j][i])
					_repr_ += '\t'
				_repr_ += '\n'
			_repr_ = _repr_.strip()
			return(_repr_)
		else:
			return(self.table.__repr__())
	def __str__(self):
		return(self.table.__str__())

class uio_struct(dict):
	__realna	= 'real'
	__intena	= 'integer'
	__charna	= 'character'
	__compna	= 'complex'
	__tabna		= 'table'
	__content	= ''
	__tabs_type	= 20
	__tabs_shape	= 10
	__pmd_digits	= 3
	def dq_units(self):
		units = dict()
		if hasattr(self, 'dq') and hasattr(self, 'z'):
			if type(self.dq)==uio_struct and type(self.z)==uio_struct:
				units = dict(
						T	= ('Temperature','K'),
						P	= ('Pressure',self.z['ei'].unit+'*'+self.z['rho'].unit),
						s	= ('Entropy',self.z['ei'].unit+'/K'),
						j1	= ('Current density, 1. component','G*'+self.z['v1'].unit),
						j2	= ('Current density, 2. component','G*'+self.z['v2'].unit),
						j3	= ('Current density, 3. component','G*'+self.z['v3'].unit),
						jabs	= ('Absolute current density','G*'+self.z['v1'].unit),
						kappa	= ('Opacity','1/'+self.z['xb1'].unit+'/('+self.z['rho'].unit+')'),
						Bc1	= ('Magnetic field, 1. component','G'),
						Bc2	= ('Magnetic field, 2. component','G'),
						Bc3	= ('Magnetic field, 3. component','G'),
						divB	= ('Magnetic field divergence','G/'+self.z['xb1'].unit),
						vabs	= ('Absolute velocity',self.z['v1'].unit),
						vh	= ('Absolute horizontal velocity',self.z['v1'].unit),
						ekin	= ('Kinetic energy',self.z['ei'].unit),
						plin	= ('Momentum',self.z['v1'].unit+'*'+self.z['rho'].unit),
						vmflux	= ('Vertical mass flux',self.z['v3'].unit+'*'+self.z['rho'].unit),
						g1	= ('First adiabatic coefficient',''),
						g3	= ('Second adiabatic coefficient',''),
						cs	= ('Sound speed',self.z['v1'].unit),
						mach	= ('Mach number',self.z['ei'].unit),
						mmu	= ('Mean molecular mass',self.z['rho'].unit+'*('+self.z['xb1'].unit+')^3'),
						tau	= ('Optical depth',''),
						Babs	= ('Absolute magnetic field','G'),
						Bh	= ('Absolute horizontal magnetic field','G'),
						B2	= ('Signed absolute magnetic field','G^2'),
						emag	= ('Magnetic energy',self.z['ei'].unit),
						cA	= ('Alfven speed',self.z['v1'].unit),
						beta	= ('Plasma beta',self.z['ei'].unit+'*'+self.z['rho'].unit+'/B^2'),
						csca	= ('Sound speed to Alfven speed ratio',''))
		return units
	def record(self, mode='', ident='', col=''):
		"""
	This function is mainly intended to be used as a call-back function
	from Fortran routines. Its first parameter is the mode in which it
	has to be used:

	  - rename:	Renames variable inside uio_struct
	  - delete:	Remove variable inside uio_struct
	  - create:	Creates uio_struct variable inside uio_struct
	  - up:		Option intended only for call-back
	  - down:	Option intended only for call-back
	  - copy:	Option intended only for call-back
	  - link:	Option intended only for call-back

	Parameters
	----------
	mode: Mode
	ident: First argument
	col: Second argument

	First and second argument depend on mode. For example 'rename'
	requires both (old_name, new_name) but delete only requires the
	first one, etc.
		"""
		if (self._loc != ''):
			if (mode=='up' and self[self._loc]._loc==''):
				self._loc=''
				return
			if (mode=='rename' and self[self._loc]._loc==''):
				col0 = ''
				if ident != '' :
					ident0=ident
					col0=col
				else :
					labels = uio_reader_module.label.tostring().strip()
					ident0 = ''.join(labels.split(':')[0])
					if labels.find(':') != -1 :
						col0 = ''.join(labels.split(':')[1])
				if col0 == '' :
					col0 = self._loc
					self[ident0] = self[col0]
					del(self[col0])
					self._loc = ident0
					return
			self[self._loc].record(mode, ident, col)
			if (mode=='top' and self[self._loc]._loc==''):
				self._loc=''
				return
			return
		if   mode == 'up' : return
		elif mode == 'create' :
			if ident != '' : ident0=ident
			else : ident0 = uio_reader_module.label.tostring().strip()
			if not hasattr(self, ident0) :
				self[ident0]=uio_struct()
				self[ident0]._parent=self
			return
		elif mode == 'delete' :
			if ident != '' :
				if not hasattr(self, ident) :
					return
			if ident == '':
				keys = [key for key in self.viewkeys()]
				for key in keys:
					if key == '_parent': continue
					if (self.__class__ ==
					   self[key].__class__) :
						self[key].record('delete')
						del(self[key])
					elif (self[key].__class__ ==
					   uio_struct_item) :
						del(self[key])
			elif self.__class__ == self[ident].__class__ :
				self[ident].record('delete')
				del(self[ident])
			elif self[ident].__class__ == uio_struct_item :
				del(self[ident])
			return
		elif mode == 'rename' :
			if ident != '' :
				ident0=ident
				col0=col
			else :
				labels = uio_reader_module.label.tostring().strip()
				ident0 = ''.join(labels.split(':')[0])
				if labels.find(':') != -1 :
					col0 = ''.join(labels.split(':')[1])
				else: col0 = ''
			if col0 != '' :
				self[col0] = self[ident0]
				del(self[ident0])
		elif mode == 'down' :
			if ident != '' : ident0=ident
			else : ident0 = uio_reader_module.label.tostring().strip()
			if not hasattr(self, ident0) :
				self[ident0]=uio_struct()
				self[ident0]._parent=self
			if self.__class__ == self[ident0].__class__ :
				self._loc=ident0
			return
		elif mode == 'copy' :
			label = uio_reader_module.label.tostring().strip()
			name = uio_reader_module.name0.tostring().strip()
			unit = uio_reader_module.unit0.tostring().strip()
			qtype = uio_reader_module.qtype.tostring().strip()
			ndim = int(uio_reader_module.ndim)
			m1 = int(uio_reader_module.m1)
			n1 = int(uio_reader_module.n1)
			m2 = int(uio_reader_module.m2)
			n2 = int(uio_reader_module.n2)
			m3 = int(uio_reader_module.m3)
			n3 = int(uio_reader_module.n3)
			m4 = int(uio_reader_module.m4)
			n4 = int(uio_reader_module.n4)
			item = uio_struct_item(self.record)
			item.name = label
			item.desc = name
			item.unit = unit
			#print('Py: '+str(label)+' '+str(qtype))
			vstrip = np.vectorize(lambda s : s.strip(), otypes=[object])
			if ndim == 1 :
			  if n1 <= m1 :
			    if qtype == self.__intena:
			        item.value = uio_reader_module.inte0.item()
			    elif qtype == self.__realna:
			        item.value = uio_reader_module.real0.item()
			    elif qtype == self.__compna:
			        item.value = uio_reader_module.comp0.item()
			    elif qtype == self.__charna:
			        item.value = ''.join(uio_reader_module.char0).strip()
			  else :
			    if qtype == self.__intena:
			        item.value = np.array(uio_reader_module.inte1D)
			    elif qtype == self.__realna:
			        item.value = np.array(uio_reader_module.real1D)
			    elif qtype == self.__compna:
			        item.value = np.array(uio_reader_module.comp1D)
			    elif qtype == self.__charna:
			        item.value = vstrip(np.apply_along_axis(''.join,-1,uio_reader_module.char1D))
			elif ndim == 2 :
			  if qtype == self.__intena:
			      item.value = np.array(uio_reader_module.inte2D)
			  elif qtype == self.__realna:
			      item.value = np.array(uio_reader_module.real2D)
			  elif qtype == self.__compna:
			      item.value = np.array(uio_reader_module.comp2D)
			  elif qtype == self.__charna:
			        item.value = vstrip(np.apply_along_axis(''.join,-1,uio_reader_module.char2D))
			elif ndim == 3 :
			  if qtype == self.__intena:
			      item.value = np.array(uio_reader_module.inte3D)
			  elif qtype == self.__realna:
			      item.value = np.array(uio_reader_module.real3D)
			  elif qtype == self.__compna:
			      item.value = np.array(uio_reader_module.comp3D)
			  elif qtype == self.__charna:
			        item.value = vstrip(np.apply_along_axis(''.join,-1,uio_reader_module.char3D))
			elif ndim == 4 :
			  if qtype == self.__intena:
			      item.value = np.array(uio_reader_module.inte4D)
			  elif qtype == self.__realna:
			      item.value = np.array(uio_reader_module.real4D)
			  elif qtype == self.__compna:
			      item.value = np.array(uio_reader_module.comp4D)
			  elif qtype == self.__charna:
			        item.value = vstrip(np.apply_along_axis(''.join,-1,uio_reader_module.char4D))
			item._parent = self
			self[label] = item
			return
		elif mode == 'link' :
			label = uio_reader_module.label.tostring().strip()
			name = uio_reader_module.name0.tostring().strip()
			unit = uio_reader_module.unit0.tostring().strip()
			qtype = uio_reader_module.qtype.tostring().strip()
			label = uio_reader_module.label.tostring().strip()
			if hasattr(self,label) :
				item = self[label]
			else :
				item = uio_struct_item(self.record)
				item.name = label
				item.desc = name
				item.unit = unit
			item._uio_struct_item__link=False
			item.value = 'linked'
			label = uio_reader_module.label.tostring().strip()
			item._uio_struct_item__link=qtype
			item._parent = self
			self[label] = item
			return
		elif mode == 'table' :
			label = uio_reader_module.label.tostring().strip()
			name = uio_reader_module.name0.tostring().strip()
			item = uio_struct_item(self.record)
			item.name = label
			item.desc = name
			item.value = uio_table()
			self[label] = item
		elif mode.startswith('table:') :
			table = ''.join(mode.split(':')[1:]).strip()
			label = uio_reader_module.label.tostring().strip()
			name = uio_reader_module.name0.tostring().strip()
			unit = uio_reader_module.unit0.tostring().strip()
			qtype = uio_reader_module.qtype.tostring().strip()
			self[table].value.headers += [label]
			self[table].value.desc += [name]
			self[table].value.units += [unit]
			vstrip = np.vectorize(lambda s : s.strip(), otypes=[object])
			if qtype == self.__intena:
				value = np.array(uio_reader_module.inte1D)
			elif qtype == self.__realna:
				value = np.array(uio_reader_module.real1D)
			elif qtype == self.__compna:
				value = np.array(uio_reader_module.comp1D)
			elif qtype == self.__charna:
			        value = vstrip(np.apply_along_axis(''.join,-1,uio_reader_module.char1D))
			else : value = np.array([], dtype=object)
			self[table].value.table += [value]
		return
	def __init__(self, *args, **kwargs):
		super(uio_struct, self).__init__(*args, **kwargs)
		#self.__dict__ = self
		object.__setattr__(self, '__dict__', self)
		_self = self
		self._loc = ''
		self.__content = 'value'
	def __set_content__(self, content):
		self.__content = content
	@property
	def names(self):
		"""
	Sets the working mode to 'names', giving short descriptions for each
	variable.

	See 'objects' mode for deeper information.
		"""
		self.__content = 'desc'
	@property
	def values(self):
		"""
	Sets the working mode to 'values' (default mode), allowing to work
	with the contents of each variable.

	See 'objects' mode for deeper information.
		"""
		self.__content = 'value'
	@property
	def units(self):
		"""
	Sets the working mode to 'units', the units of each variable.

	See 'objects' mode for deeper information.
		"""
		self.__content = 'unit'
	@property
	def objects(self):
		"""
	Sets the working mode to 'objects'. All variables contained in
	uio_struct variables are in fact of type uio_struct_item, and the
	uio_struct_item is a class containing the value of the variable, its
	units, a short descriptions...
	
	This mode provides transparency so that variables appear as they
	really are. For instance in the (default) 'values' mode, entries in
	uio_struct variables are in fact 'links' to the value of the
	corresponding uio_struct_item variable.

	Example:

	>>> model = uio_struct()
	>>> model.theanswer = 42
	>>> type(model.theanswer)
	<type 'int'>
	>>> type(model['theanswer'])
	<class 'uio_struct_item'>
	>>> model.objects
	>>> type(model.theanswer)
	>>> <class 'uio_struct_item'>
		"""
		self.__content = 'self'
	def __setattr__(self, name, value):
		if value.__class__ != self.__class__ and value.__class__ != uio_struct_item and name[0] != '_':
			if hasattr(self, name):
				self[name].__setattr__(self.__content, value)
			else:
				item = uio_struct_item(self.record)
				item.__setattr__(self.__content, value)
				self[name] = item
				self.__dict__[name] = item
		else :
			self[name] = value
			if not name.startswith('_'):
				self.__dict__[name] = value

	def __getattribute__(self, name):
		attr = object.__getattribute__(self, name)
		if attr.__class__ == uio_struct_item :
			#return attr
			#return object.__getattribute__(attr, 'value')
			return attr.__getattribute__(self.__content)
		elif attr.__class__ == uio_struct :
			attr.__set_content__(self.__content)
			return attr
		else :
			return attr
	def __repr__(self):
		_repr_ = ''
		for key in sorted(self.viewkeys()):
			if str(key)[0] == '_': continue
			if type(self[key]) == uio_struct_item :
			#if hasattr(object.__getattribute__(self,key), 'value'):
				link = self[key]._uio_struct_item__link
				self[key]._uio_struct_item__link = False
				if type(self[key].value) == str:
					if self[key].value == 'linked':
						self[key]._uio_struct_item__link = link
				_repr_ += (str(key)+'\t'+str(self[key].value.__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key].value, 'shape'):
					if self[key].value.shape != ():
						_repr_ += '\t'+str(self[key].value.shape)
						_repr_ = _repr_.expandtabs(self.__tabs_shape)
				self[key]._uio_struct_item__link = link
			else:
				_repr_ += (str(key)+'\t'+str(self[key].__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key], 'shape'):
					_repr_ += '\t'+str(self[key].shape)
					_repr_ = _repr_.expandtabs(self.__tabs_shape)
			#_repr_ += str(key)+'\t'+str(self[key].__class__)
			#if hasattr(self[key], 'shape'):
				#_repr_ += '\t'+str(self[key].shape)
			_repr_ += '\n'
		_repr_ = _repr_.strip()
		return(_repr_)
	def __str__(self):
		return(super(uio_struct, self).__repr__())
	def __map__(self, depth=float('inf')):
		_map_ = ''
		for key in sorted(self.viewkeys()):
			if str(key)[0] == '_': continue
			if type(self[key]) == uio_struct_item :
			#if hasattr(object.__getattribute__(self,key), 'value'):
				link = self[key]._uio_struct_item__link
				self[key]._uio_struct_item__link = False
				if type(self[key].value) == str:
					if self[key].value == 'linked':
						self[key]._uio_struct_item__link = link
				_map_ += (str(key)+'\t'+str(self[key].value.__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key].value, 'shape'):
					if self[key].value.shape != ():
						_map_ += '\t'+str(self[key].value.shape)
						_map_ = _map_.expandtabs(self.__tabs_shape)
				self[key]._uio_struct_item__link = link
			else:
				_map_ += (str(key)+'\t'+str(self[key].__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key], 'shape'):
					_map_ += '\t'+str(self[key].shape)
					_map_ = _map_.expandtabs(self.__tabs_shape)
			_map_ += '\n'
			if depth>0:
				#_val_ = str(self[key].__repr__())
				_val_ = repr(self.__getattribute__(key))
				if self.__class__ == self[key].__class__ :
					_val_ = self[key].__map__(depth-1)
					_val_ = _val_.replace('\n','\n\t')
					_val_ = '\t'+_val_
				_map_ += _val_
				_map_ += '\n\n'
		_map_ = _map_.strip()
		return(_map_)
	def map(self, depth=float('inf')):
		"""
	Prints a map of an uio_struct variable. The map is by default
	recursive; one can limit its depth by setting the depth argument
	to some finite value.

	Note that the output depends on the mode:

	 - values:	displays variables values
	 - units:	displays variables units
	 - names:	displays variables names (short description)
	 - objects:	displays variables values
		"""
		_map_ = ''
		for key in sorted(self.viewkeys()):
			if str(key)[0] == '_': continue
			if type(self[key]) == uio_struct_item :
			#if hasattr(object.__getattribute__(self,key), 'value'):
				link = self[key]._uio_struct_item__link
				self[key]._uio_struct_item__link = False
				if type(self[key].value) == str:
					if self[key].value == 'linked':
						self[key]._uio_struct_item__link = link
				_map_ += (str(key)+'\t'+str(self[key].value.__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key].value, 'shape'):
					if self[key].value.shape != ():
						_map_ += '\t'+str(self[key].value.shape)
						_map_ = _map_.expandtabs(self.__tabs_shape)
				self[key]._uio_struct_item__link = link
			else:
				_map_ += (str(key)+'\t'+str(self[key].__class__)).expandtabs(self.__tabs_type)
				if hasattr(self[key], 'shape'):
					_map_ += '\t'+str(self[key].shape)
					_map_ = _map_.expandtabs(self.__tabs_shape)
				self[key]._uio_struct_item__link = link
			#_map_ += str(key)+'\t'+str(self[key].__class__)
			#if hasattr(self[key], 'shape'):
			#	_map_ += '\t'+str(self[key].shape)
			_map_ += '\n'
			if depth>0:
				#_val_ = str(self[key].__repr__())
				_val_ = repr(self.__getattribute__(key))
				if self.__class__ == self[key].__class__ :
					_val_ = self[key].__map__(depth-1)
					_val_ = _val_.replace('\n','\n\t')
					_val_ = '\t'+_val_
				_map_ += _val_
				_map_ += '\n\n'
		_map_ = _map_.strip()
		print(_map_)
		return
	def load(self, *args, **kwargs):
		"""
	Loads UIO data into uio_struct variable. Opens model file, reads data,
	stores it an closes model file again.
	
	If both fullfile and parfile are provided, EOS and OPA data are also
	loaded from data files specified in parfile in order to provided
	derived quantities.

	Parameters
	----------
	modelfile : Model file: FULL/MEAN/STA/END file.
	parfile : Parameter file (only when fullfile is provided).
	imodel : Model index (default 0) when model file contains multiple
	         snapshots.
		"""
		self.record('top')
		modelfile = ''
		parfile = ''
		imodel = None
		if len(args) >= 1 :
			modelfile = args[0]
		if len(args) == 2 :
			if type(args[1]) == str :
				parfile = args[1]
			elif type(args[1]) == int :
				imodel = args[1]
			else :
				raise ValueError('Second argument must be of type int or str.')
		elif len(args) >= 3 :
			parfile = args[1]
			imodel = args[2]
		for key in kwargs :
			if key == 'modelfile' and modelfile == '' :
				modelfile = kwargs[key]
			elif key == 'parfile' and parfile == '' :
				parfile = kwargs[key]
			elif key == 'n' and imodel == None :
				imodel = kwargs[key]
		if imodel == 0 or imodel == None : imodel = 1
		if with_pypmd:
			pmd_handle = pypmd.openrd(modelfile)
			if pmd_handle[-1] == -500 :
				pypmd.close(pmd_handle[0])
				pyread(self.record, modelfile, parfile, imodel)
				self._datatype = 'uio'
			elif pmd_handle[-1] >= 0 :
				self.pmd_parse(pmd_handle)
				self._datatype = pmd_handle[1]
			else :
				pypmd.close(pmd_handle[0])
		else:
			pyread(self.record, modelfile, parfile, imodel)
			self._datatype = 'uio'
		if parfile != '' :
			if not hasattr(self, 'dq') :
				self.record('create', 'dq')
				self.dq.box_id = 'dq'
			dqs = self.dq_units()
			for dq in dqs:
				self.dq.__setattr__(dq, None)
				self.dq[dq].name  = dq
				self.dq[dq].desc  = dqs[dq][0]
				self.dq[dq].unit  = dqs[dq][1]
				self.dq[dq]._uio_struct_item__link='real'
	def open(self, *args, **kwargs):
		"""
	Opens model file for reading. This function is only used in
	combination with header, next, count, imodel, and close when the
	number of snapshots in the model file is unknown or when one
	wishes to read all snapshots.

	In any other circumstance use the load function instead.

	If both fullfile and parfile are provided, EOS and OPA data are also
	loaded from data files specified in parfile in order to provided
	derived quantities.

	Parameters
	----------
	modelfile : Model file: FULL/MEAN/STA/END file.
	parfile : Parameter file (only when fullfile is provided).
		"""
		self.record('top')
		modelfile = ''
		parfile = ''
		if len(args) >= 1 :
			modelfile = args[0]
		if len(args) == 2 :
			parfile = args[1]
		for key in kwargs :
			if key == 'modelfile' and modelfile == '' :
				modelfile = kwargs[key]
			elif key == 'parfile' and parfile == '' :
				parfile = kwargs[key]
		open_model(modelfile, parfile)
	@property
	def header(self):
		"""
	Read current snapshot header in opened modelfile. Save header in
	uio_struct variable. After opening model file with open, header
	must be either read before using next, or skipped until file is
	closed and reopened.
		"""
		read_header(self.record)
	@property
	def count(self):
		"""
	Returns the number of snapshots in opened model file.
		"""
		count_model()
		return uio_reader_module.r_nmodel.item()
	@property
	def imodel(self):
		"""
	Returns the current snapshot index in opened model file.
		"""
		return uio_reader_module.r_imodel.item()
	@property
	def next(self):
		"""
	Loads next snapshot of opened model file into uio_struct variable.
	After this function is called the first time, header cannot be
	called anymore.

	Note: open will not load any data into uio_struct variable; the
	      first call to next will load the first snapshot to
	      uio_struct variable.
		"""
		keys = [key for key in self.viewkeys()]
		for key in keys:
			if key == '_parent': continue
			if (self.__class__ ==
			   self[key].__class__) :
				if key != 'head' and key != 'dq':
					self[key].record('delete')
					del(self[key])
			elif (self[key].__class__ ==
			   uio_struct_item) :
				del(self[key])
		read_next(self.record)
		if hasattr(self, 'dq') and uio_reader_module.r_imodel.item() >= 0 :
			self.record('delete','dq')
		if ''.join(uio_reader_module.r_parfile).strip() != '' and uio_reader_module.r_imodel.item() >= 0 :
			self.record('create', 'dq')
			self.dq.box_id = 'dq'
			dqs = self.dq_units()
			for dq in dqs:
				self.dq.__setattr__(dq, None)
				self.dq[dq].name  = dq
				self.dq[dq].desc  = dqs[dq][0]
				self.dq[dq].unit  = dqs[dq][1]
				self.dq[dq]._uio_struct_item__link='real'
		return uio_reader_module.r_imodel.item()
	@property
	def close(self):
		"""
	Closes model file. This function is only intended to manually close
	a file opened with open.

	Note: after reading a model with load, the file is automatically
	closed.
		"""
		close_model()
	def export(self, filename, **kwargs):
	        # Parse optional arguments
        	for key in kwargs:
	                if key=='max_sz': max_sz=kwargs[key]
	                elif key=='box': box=kwargs[key]
	                elif key=='centre': centre=kwargs[key]
	                elif key=='dtype': arrtype=kwargs[key]
	        if not 'max_sz' in locals(): max_sz=(-1,-1,-1)
	        if not 'box' in locals(): box=(-1,-1,-1)
	        if not 'centre' in locals(): centre=(-1,-1,-1)
	        if not 'dtype' in locals(): arrtype=np.float32

		nX = np.shape(self.z['xb1'].value)[0]
		nY = np.shape(self.z['xb2'].value)[1]
		nZ = np.shape(self.z['xb3'].value)[2]
		periodX = self.z['xb1'].value[-1,0,0]-self.z['xb1'].value[0,0,0]
		periodY = self.z['xb2'].value[0,-1,0]-self.z['xb2'].value[0,0,0]
                id_min0=np.array([0,0,0])
                id_max0=np.array([nX, nY, nZ])
                id_min=id_min0.copy()
                id_max=id_max0.copy()
                steps=np.array([1,1,1])
		periodicity_flag = False
		periodicity_dir = [0,0]
                for i in range(3):
                        if centre[i]<0 and box[i]>=0 : centre[i]=0
                        if centre[i]>=0 :
                                if box[i]>=0 :
                                        #id_min[i]=max(id_min0[i],int(np.ceil(centre[i]-box[i]/2.)))
                                        #id_max[i]=min(id_max0[i],1+int(np.floor(centre[i]+box[i]/2.)))
					if i < 2:
	                                        id_min[i]=int(np.ceil(centre[i]-box[i]/2.))-1
	                                        id_max[i]=int(np.ceil(centre[i]+box[i]/2.))
						if ( id_min0[i] > int(np.ceil(centre[i]-box[i]/2.))-1 ):
							periodicity_flag = True
							periodicity_dir[i] = -1
						if ( id_max0[i] < int(np.ceil(centre[i]+box[i]/2.)) ):
							periodicity_flag = True
							periodicity_dir[i] = 1
					else:
	                                        id_min[i]=max(id_min0[i],int(np.ceil(centre[i]-box[i]/2.))-1)
	                                        id_max[i]=min(id_max0[i],int(np.ceil(centre[i]+box[i]/2.)))
                                        if max_sz[i]>=0 :
                                                steps[i]=max(1,int(np.floor((id_max[i]-id_min[i]+1)/max_sz[i])))
                                elif max_sz[i]>=0 :
                                        #id_min[i]=max(id_min[i],int(np.ceil(centre[i]-max_sz[i]/2.)))
                                        #id_max[i]=min(id_max[i],1+int(np.floor(centre[i]+max_sz[i]/2.)))
					if i < 2:
                                                id_min[i]=int(np.ceil(centre[i]-max_sz[i]/2.))-1
                                                id_max[i]=int(np.ceil(centre[i]+max_sz[i]/2.))
						if ( id_min0[i] > int(np.ceil(centre[i]-max_sz[i]/2.))-1 ):
							periodicity_flag = True
							periodicity_dir[i] = -1
						if ( id_max0[i] < int(np.ceil(centre[i]+max_sz[i]/2.)) ):
							periodicity_flag = True
							periodicity_dir[i] = 1
					else:
                                                id_min[i]=max(id_min0[i],int(np.ceil(centre[i]-max_sz[i]/2.))-1)
                                                id_max[i]=min(id_max0[i],int(np.ceil(centre[i]+max_sz[i]/2.)))
                        else :
                                if max_sz[i]>=0 :
                                        steps[i]=max(1,int(np.floor((id_max[i]-id_min[i]+1)/max_sz[i])))
		m1=id_min[0]+1
		m2=id_min[1]+1
		m3=id_min[2]+1
		n1=id_max[0]-1
		n2=id_max[1]-1
		n3=id_max[2]-1
		meshX, meshY, meshZ = np.mgrid[id_min[0]:id_max[0]:steps[0],id_min[1]:id_max[1]:steps[1],id_min[2]:id_max[2]:steps[2]]
		#meshX = meshX.T
		#meshY = meshY.T
		#meshZ = meshZ.T
		if hasattr(self.z, 'bb1'): Bb_flag = True
		else : Bb_flag = False
		if steps[0]==1 and not periodicity_flag :
			xb1=np.asfortranarray(self.z['xb1'].value[id_min[0]:id_max[0],:,:])
			xc1=np.asfortranarray(self.z['xc1'].value[id_min[0]:id_max[0]-1,:,:])
                else :
                        #xb1=np.array(self.z['xb1'].value[id_min[0]:id_max[0]:steps[0],:,:],order='F', dtype=arrtype)
                        #xc1=np.array(self.z['xc1'].value[id_min[0]:id_max[0]-1:steps[0],:,:],order='F', dtype=arrtype)
			rangeX = np.arange(id_min[0],id_max[0],steps[0])%(nX-1)
			xb1=np.array(self.z['xb1'].value[rangeX,:,:],order='F', dtype=arrtype)
			xc1=np.array(self.z['xc1'].value[rangeX[:-1],:,:],order='F', dtype=arrtype)
			xb1 = extend_periodically(xb1,periodX,periodicity_dir[0])
			xc1 = extend_periodically(xc1,periodX,periodicity_dir[0])
                if steps[1]==1 and not periodicity_flag :
			xb2=np.asfortranarray(self.z['xb2'].value[:,id_min[1]:id_max[1],:])
			xc2=np.asfortranarray(self.z['xc2'].value[:,id_min[1]:id_max[1]-1,:])
                else :  
                        #xb2=np.array(self.z['xb2'].value[:,id_min[1]:id_max[1]:steps[1],:],order='F', dtype=arrtype)
                        #xc2=np.array(self.z['xc2'].value[:,id_min[1]:id_max[1]-1:steps[1],:],order='F', dtype=arrtype)
			rangeY = np.arange(id_min[1],id_max[1],steps[1])%(nY-1)
			xb2=np.array(self.z['xb2'].value[:,rangeY,:],order='F', dtype=arrtype)
			xc2=np.array(self.z['xc2'].value[:,rangeY[:-1],:],order='F', dtype=arrtype)
			xb2 = extend_periodically(xb2,periodY,periodicity_dir[1])
			xc2 = extend_periodically(xc2,periodY,periodicity_dir[1])
                if steps[2]==1 and not periodicity_flag :
			xb3=np.asfortranarray(self.z['xb3'].value[:,:,id_min[2]:id_max[2]])
			xc3=np.asfortranarray(self.z['xc3'].value[:,:,id_min[2]:id_max[2]-1])
                else :
                        #xb3=np.array(self.z['xb3'].value[:,:,id_min[2]:id_max[2]:steps[2]],order='F', dtype=arrtype)
                        #xb3=np.array(self.z['xc3'].value[:,:,id_min[2]:id_max[2]-1:steps[2]],order='F', dtype=arrtype)
			rangeZ = np.arange(id_min[2],id_max[2],steps[2])%(nZ-1)
			xb3=np.array(self.z['xb3'].value[:,:,rangeZ],order='F', dtype=arrtype)
			xc3=np.array(self.z['xc3'].value[:,:,rangeZ[:-1]],order='F', dtype=arrtype)
		v1=self.z['v1'].value
		v2=self.z['v2'].value
		v3=self.z['v3'].value
		rho=self.z['rho'].value
		ei=self.z['ei'].value
		boxes = ['v1', 'v2', 'v3', 'rho', 'ei'];
		if Bb_flag:
			Bb1=self.z['bb1'].value
			Bb2=self.z['bb2'].value
			Bb3=self.z['bb3'].value
			B1_unit=self.z['bb1'].unit
			B2_unit=self.z['bb2'].unit
			B3_unit=self.z['bb3'].unit
		else :
			Bb1=np.zeros((nX,nY-1,nZ-1),dtype=arrtype,order='F')
			Bb2=np.zeros((nX-1,nY,nZ-1),dtype=arrtype,order='F')
			Bb3=np.zeros((nX-1,nY-1,nZ),dtype=arrtype,order='F')
			B1_unit=''
			B2_unit=''
			B3_unit=''
                if (not np.all(steps==[1,1,1])) :
			m1=1
			m2=1
			m3=1
			n1=id_max[0]-id_min[0]-1
			n2=id_max[1]-id_min[1]-1
			n3=id_max[2]-id_min[2]-1
                #if (not np.all(steps==[1,1,1])) or periodicity_flag :
		for box in boxes:
			#vars()[box] = vars()[box][id_min[0]:id_max[0]-1:mesh_step[mesh_id0[index]][0],id_min[1]:id_max[1]-1:mesh_step[mesh_id0[index]][1],id_min[2]:id_max[2]-1:mesh_step[mesh_id0[index]][2]]
			quantitiy = np.asfortranarray(vars()[box][(meshX%(nX-1))[:-1,:-1,:-1],(meshY%(nY-1))[:-1,:-1,:-1],(meshZ%(nZ-1))[:-1,:-1,:-1]])
			exec(box+'=quantitiy')
		Bb1 = np.asfortranarray(Bb1[(meshX%(nX-1))[:,:-1,:-1],(meshY%(nY-1))[:,:-1,:-1],(meshZ%(nZ-1))[:,:-1,:-1]])
		Bb2 = np.asfortranarray(Bb2[(meshX%(nX-1))[:-1,:,:-1],(meshY%(nY-1))[:-1,:,:-1],(meshZ%(nZ-1))[:-1,:,:-1]])
		Bb3 = np.asfortranarray(Bb3[(meshX%(nX-1))[:-1,:-1,:],(meshY%(nY-1))[:-1,:-1,:],(meshZ%(nZ-1))[:-1,:-1,:]])
		if hasattr(self.z,'time_db'):
			time_db=self.z['time_db'].value
		else: time_db=self['modeltime'].value
		if hasattr(self, 'head'):
			description=self.head['description'].value
			version=self.head['version'].value
			history=self.head['history'].value
		else :
			description=''
			version='Unknown version'
			history=np.array(['Unknown past history'])
		if hasattr(self, 'time_out_mean_last') :
			toml = self['time_out_mean_last'].value
		else :
			toml = 0.0
		if hasattr(self, 'time_out_full_last') :
			tofl = self['time_out_full_last'].value
		else :
			tofl = 0.0
		history=np.array([list(h.ljust(80)) for h in history])
		#print('Filename: '+filename)
		#print('xb1 shape: '+str(np.shape(xb1)))
		#print('xb2 shape: '+str(np.shape(xb2)))
		#print('xb3 shape: '+str(np.shape(xb3)))
		#print('xc1 shape: '+str(np.shape(xc1)))
		#print('xc2 shape: '+str(np.shape(xc2)))
		#print('xc3 shape: '+str(np.shape(xc3)))
		#print('v1 shape: '+str(np.shape(v1)))
		#print('v2 shape: '+str(np.shape(v2)))
		#print('v3 shape: '+str(np.shape(v3)))
		#print('rho shape: '+str(np.shape(rho)))
		#print('ei shape: '+str(np.shape(ei)))
		#print('Bb_flag: '+str(Bb_flag))
		#print('Bb1 shape: '+str(np.shape(Bb1)))
		#print('Bb2 shape: '+str(np.shape(Bb2)))
		#print('Bb3 shape: '+str(np.shape(Bb3)))
		#print('B1_unit: '+B1_unit)
		#print('B2_unit: '+B2_unit)
		#print('B3_unit: '+B3_unit)
		#print('dtime: '+str(self['dtime'].value))
		#print('modelitime: '+str(self['modelitime'].value))
		#print('modeltime: '+str(self['modeltime'].value))
		#print('time_db: '+str(time_db))
		#print('history: '+str(history))
		#print('version: '+version)
		#print('toml: '+str(toml))
		#print('tofl: '+str(tofl))
		#print('m1: '+str(m1))
		#print('n1: '+str(n1))
		#print('m2: '+str(m2))
		#print('n2: '+str(n2))
		#print('m3: '+str(m3))
		#print('n3: '+str(n3))
		write_model(filename, xb1, xb2, xb3, xc1, xc2, xc3, v1, v2, v3, rho, ei, Bb_flag, Bb1, Bb2, Bb3, B1_unit, B2_unit, B3_unit, self['dtime'].value, self['modelitime'].value, self['modeltime'].value, time_db, description, history, version, toml, tofl, m1, n1, m2, n2, m3, n3)
	# PMD routines
	def pmd_parse(self, pmd_handle):
		unit			= pmd_handle[0]
		self.magic_str		= pmd_handle[1]
		self.endianness		= pmd_handle[2]
		self.int_sz		= pmd_handle[3]
		self.db_sz		= pmd_handle[4]
		self.pmd_version	= pmd_handle[5]
		self.localtime		= pmd_handle[6]
		self.periodic		= pmd_handle[7]
		self.domain_sz		= pmd_handle[8]
		self.domain_origin	= pmd_handle[9]
		self.dimensions		= pmd_handle[10]
		self.xc1		= pmd_handle[11]
		self.xc2		= pmd_handle[12]
		self.xc3		= pmd_handle[13]
		self.n_radtheta		= pmd_handle[14]
		self.n_radphi		= pmd_handle[15]
		self.atomic_module	= pmd_handle[16].strip()
		self.pmd_comment	= pmd_handle[17].strip()
		self.module_hd_sz	= pmd_handle[18]
		self.unused_var		= pmd_handle[19]
		#ndatasets		= pmd_handle[20]
		node_sz			= pmd_handle[20]
		if with_pypmd:
			module_header = list(pypmd.rd_module_hd(unit,self.module_hd_sz))
			for i in range(len(module_header)):
				if module_header[i] == '':
					module_header[i] = '\x00'
			self.module_header	= bytes(''.join(module_header))
			#for i in range(1,ndatasets+1):
			#	data='data'+('{0:0>'+str(self.__pmd_digits)+'}').format(i)
			#	self[data] = uio_struct_item()
			#	self[data].value = np.empty(self.dimensions, order='F')
			#	pypmd.rd_box(unit, self[data].value)
			N = np.prod(self.dimensions)
			nx, ny, nz = self.dimensions
			node_data = []
			for i in range(N):
				data = list(pypmd.rd_node(unit, node_sz))
				for j in range(len(data)):
					if data[j] == '':
						data[j] = '\x00'
					data[j] = ord(data[j])
				#node_data += [bytes(''.join(data))]
				node_data += [data]
			print(np.size(node_data))
			print(str(node_sz))
			self.node_data = np.transpose(np.array(node_data, dtype=np.uint8).reshape(nz,ny,nx,node_sz),(2,1,0,3))
			pypmd.close(unit)
	def pmd_write(self, filename):
		if with_pypmd:
			xc1 = np.empty(8192)
			xc2 = np.empty(8192)
			xc3 = np.empty(8192)
			xc1[0:len(self.xc1)] = self.xc1
			xc2[0:len(self.xc2)] = self.xc2
			xc3[0:len(self.xc3)] = self.xc3
			self.xc1 = xc1
			self.xc2 = xc2
			self.xc3 = xc3
			unit, self.localtime, stat = pypmd.openwr(filename, self.magic_str, self.endianness, self.int_sz, self.db_sz, self.pmd_version, self.localtime, self.periodic, self.domain_sz, self.domain_origin, self.dimensions, self.xc1, self.xc2, self.xc3, self.n_radtheta, self.n_radphi, self.atomic_module, self.pmd_comment, self.module_hd_sz, self.unused_var)
			if self.module_hd_sz > 0:
				if len(self.module_header)<self.module_hd_sz:
					self.module_header = self.module_header+(self.module_hd_sz-len(self.module_header))*'\x00'
				self.module_header = self.module_header[0:self.module_hd_sz]
				pypmd.wr_module_hd(unit, list(self.module_header))
			for i in range(1,10**self.__pmd_digits-1):
				data='data'+('{0:0>'+str(self.__pmd_digits)+'}').format(i)
				if hasattr(self, data):
					pypmd.append_box(unit, self[data].value)
				else:
					break
			if hasattr(self, 'node_data'):
				N = np.prod(self.dimensions)
				node_sz = int(np.size(self.node_data)/N)
				node_data = np.transpose(self.node_data,(2,1,0,3)).reshape((N,node_sz))
				for node in node_data:
					node = [chr(i) for i in node]
					pypmd.append_node(unit, node)
			pypmd.close(unit)
		else:
			print("You need to recompile with pypmd support.")

def pmd_struct():
	pmdstruct		= uio_struct()
	pmdstruct.magic_str	= "portapmd"
	pmdstruct.endianness	= 0
	pmdstruct.int_sz	= 4
	pmdstruct.db_sz		= 8
	pmdstruct.pmd_version	= 2
	pmdstruct.localtime	= np.array([0, 0, 0, 0, 0, -1])
	pmdstruct.periodic	= np.array([0, 0])
	pmdstruct.domain_sz	= np.array([0., 0., 0.])
	pmdstruct.domain_origin	= np.array([0., 0., 0.])
	pmdstruct.dimensions	= np.array([0, 0, 0])
	pmdstruct.xc1		= np.zeros(8192)
	pmdstruct.xc2		= np.zeros(8192)
	pmdstruct.xc3		= np.zeros(8192)
	pmdstruct.n_radtheta	= 0
	pmdstruct.n_radphi	= 0
	pmdstruct.atomic_module	= ''
	pmdstruct.pmd_comment	= ''
	pmdstruct.module_hd_sz	= 0
	pmdstruct.unused_var	= 0
	pmdstruct.module_header	= ''
	return pmdstruct
init()
