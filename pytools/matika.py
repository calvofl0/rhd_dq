#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import gradient, shape, asarray, meshgrid, arange, ones_like, zeros, compress
from numpy import sum as asum
from numpy.linalg import norm
try:
	from matplotlib.pyplot import quiver
	matplotlib_flag = True
except ImportError:
	matplotlib_flag = False
	from warnings import warn
	warn("Matplotlib is not installed", ImportWarning)

'''
This module provides general functions for tensor calculus.

The index convention is that vector fields components are ordered with
the same order as coordinates.
'''

def advection(v, f, *varargs, **kwargs):
	'''
Return advection of an N-dimensional vector field f (or a scalar field f)
along and other N-dimensional vector field v.

v is an array-like (N+1)-dimensional object containing a vector field
(the first index specifies the component of the vector field).
f is either a vector field with same shape as v, or a scalar field (ie
an N-dimensional array of same shape as each component of v).

Uses the numpy.gradient function to compute partial derivatives and
inherits the same optional arguments.
	'''
	if len(shape(v)) == 0:
		raise TypeError('Invalid vector field v')
	N = shape(v)[0]
	shp = shape(v[0])
	if len(shp) != N:
		raise TypeError('Invalid vector field v')
	D = []
	for d in v:
		if shape(d) != shp:
			raise TypeError('Invalid vector field v')
	if len(shape(f)) == N:
		D = gradient(f, *varargs, **kwargs)
		return asum(asarray(v)*D, axis=0)
	elif len(shape(f)) == N+1:
		for d in f:
			D += [gradient(d, *varargs, **kwargs)]
	else:
			raise TypeError('Invalid scalar/vector field f')

	return [asum(asarray(v)*D[i], axis=0) for i in range(N)]

def lie(f, v, *varargs, **kwargs):
	'''
Return Lie derivative of an N-dimensional vector field f (or a scalar
field f) along and other N-dimensional vector field v.

v is an array-like (N+1)-dimensional object containing a vector field
(the first index specifies the component of the vector field).
f is either a vector field with same shape as v, or a scalar field (ie
an N-dimensional array of same shape as each component of v).

Uses the numpy.gradient function to compute partial derivatives and
inherits the same optional arguments.
	'''
	if len(shape(v)) == 0:
		raise TypeError('Invalid vector field v')
	N = shape(v)[0]
	adv = advection(f, v, *varargs, **kwargs)
	if len(shape(f)) == N:
		Dlie = adv
	elif len(shape(f)) == N+1:
		adv_vf = advection(v, f, *varargs, **kwargs)
		Dlie = [adv[i] - adv_vf[i] for i in range(N)]
	else:
		raise TypeError('Invalid scalar/vector field f')

	return Dlie

def rotational(v, *varargs, **kwargs):
	'''
Return rotational of a 3-dimensional vector field.

v is an array-like (3+1)-dimensional object containing a 3D vector field
(the first index specifies the component of the vector field).

Uses the numpy.gradient function to compute partial derivatives and
inherits the same optional arguments.

Note: if the coordinates are right-handed and ordered in increasing order,
      the result is the ordinary (right-handed) rotational.
	'''
	if len(shape(v)) == 0:
		raise TypeError('Invalid vector field v')
	N = shape(v)[0]
	if N != 3:
		raise TypeError('Vector field v is not a 3D vector field')
	if len(shape(v)) != N+1:
		raise TypeError('Invalid vector field v')
	shp = shape(v[0])
	if len(shp) != N:
		raise TypeError('Invalid vector field v')
	D = []
	for d in v:
		if shape(d) != shp:
			raise TypeError('Invalid vector field v')
		D += [gradient(d, *varargs, **kwargs)]

	return [D[(i-1)%N][(i+1)%N]-D[(i+1)%N][(i-1)%N] for i in range(N)]

def mgrid(*varargs):
	'''
Return a meshgrid with the indexing convention of this module (vector
components and tensor field coordinates are equally ordered).

Coordinates along each dimension must be provided as 1-arrays or as
scalars (in which case coordinates are given by np.arange of the scalar
values).
	'''
	if len(varargs) == 0:
		raise TypeError('You must provide coordinates')
	shp = shape(varargs[0])
	#for arg in varargs:
	#	if shape(arg) != shp:
	#		raise TypeError('Invalid coordinates')
	if len(shp) == 0:
		coords = []
		for arg in varargs:
			coords += [arange(arg)]
		return meshgrid(*coords, indexing='ij')
	elif len(shp) == 1:
		return meshgrid(*varargs, indexing='ij')
	else:
		return varargs

class VectorField(object):
	'''
Implements a vector field with standard term-by-term operations:
   * Addition
   * Substraction
   * Multiplication
   * Division
Also implements other operations:
   * Cross product (for 3-D fields only, with the __xor__ function (^))

A vector is considered as a vector field with constant components.
A scalar field is considered as a vector field with the same components
along each direction
A scalar is considered as a vector field with the same constant components
along each direction
	'''
	def __init__(self, v):
		if type(v) == type(self):
			self._v = v._v
		else:
			if len(shape(v)) == 0:
				raise TypeError('Invalid vector field')
			self._v = asarray(v)
	def __add__(self, w):
		if type(w) is type(self):
			w0 = w._v
		else:
			w0 = w
		shp = shape(w0)
		if len(shp) == 0:
			return VectorField([u+w0 for u in self._v])
		elif shp == shape(self._v[0]):
			return VectorField(self._v+w0)
		elif len(self._v) == len(w0):
			return VectorField([self._v[i]+w0[i] for i in range(len(w0))])
		else:
			raise TypeError('Invalid scalar/vector')
	def __sub__(self, w):
		if type(w) is type(self):
			w0 = w._v
		else:
			w0 = w
		shp = shape(w0)
		if len(shp) == 0:
			return VectorField([u-w0 for u in self._v])
		elif shp == shape(self._v[0]):
			return VectorField(self._v-w0)
		elif len(self._v) == len(w0):
			return VectorField([self._v[i]-w0[i] for i in range(len(w0))])
		else:
			raise TypeError('Invalid scalar/vector')
	def __mul__(self, w):
		if type(w) is type(self):
			w0 = w._v
		else:
			w0 = w
		shp = shape(w0)
		if len(shp) == 0:
			return VectorField([u*w0 for u in self._v])
		elif shp == shape(self._v[0]):
			return VectorField(self._v*w0)
		elif len(self._v) == len(w0):
			return VectorField([self._v[i]*w0[i] for i in range(len(w0))])
		else:
			raise TypeError('Invalid scalar/vector')
	def __div__(self, w):
		if type(w) is type(self):
			w0 = w._v
		else:
			w0 = w
		shp = shape(w0)
		if len(shp) == 0:
			return VectorField([u/w0 for u in self._v])
		elif shp == shape(self._v[0]):
			return VectorField(self._v/w0)
		elif len(self._v) == len(w0):
			return VectorField([self._v[i]/w0[i] for i in range(len(w0))])
		else:
			raise TypeError('Invalid scalar/vector')
	def __xor__(self, w):
		if type(w) is type(self):
			w0 = w._v
		else:
			w0 = w
		shp = shape(w0)
		if len(shp) == 0 or shp == shape(self._v[0]):
			raise TypeError('Invalid vector (scalar provided)')
		elif len(self._v) == len(w0) == 3:
			v0=self._v
			return VectorField([v0[(i+1)%3]*w0[(i-1)%3]-v0[(i-1)%3]*w0[(i+1)%3] for i in range(3)])
		else:
			raise TypeError('Invalid 3-vector')
	def __iter__(self):
		return self._v.__iter__()
	def __getitem__(self, index):
		return self._v.__getitem__(index)
	def __len__(self):
		return self._v.__len__()
	def __repr__(self):
		return self._v.__repr__()
	def __str__(self):
		return self._v.__str__()
	@property
	def shape(self):
		return self._v.shape
	def quiver(self, *varargs, **kwargs):
		'''
Plots the projection of the vector field on a given plane.

Parameters
----------
plane:   integer, (0, 1 or 2, correspond to the perpendicular axis of the
                   to be plotted)
index:   integer, (index of slice to plot)
X, Y, Z: 1D arrays or meshgrids, optional (default spacing is 1)

C:       array, used to map colors to the arrows (see matplotlib doc)
kwargs:  
	 * density: scalar or list of two elements, distance between each
	            arrow in the two-dimensional plot
	 * see matplotlib doc for other options
		'''
		if not matplotlib_flag:
			raise ImportError('Matplotlib is not installed')
		if len(varargs) not in set([2,3,5,6]):
			raise TypeError('Invalid number of arguments')
		if len(self._v) != 3:
			raise TypeError('Invalid 3D vector field')
		plane = varargs[0]
		index = varargs[1]
		C = None
		coords = None
		if len(varargs) == 3:
			C = varargs[2]
		elif len(varargs) >= 5:
			coords = mgrid(*varargs[2:5])
		if len(varargs) == 6:
			C = varargs[5]
		dens = kwargs.pop('density', 1)
		if len(shape(dens)) == 0:
			dens = [dens, dens]
		elif len(shape(dens)) != 1 or len(dens) != 2:
			raise TypeError('Invalid density argument')
		shp = shape(self._v[0])
		nshp = tuple(asarray(shp)[[(plane+1)%3, (plane-1)%3]])
		cond = zeros(shp[plane], dtype=bool)
		cond[index] = True
		U = compress(cond, self._v[(plane+1)%3], axis=plane).reshape(nshp)
		V = compress(cond, self._v[(plane-1)%3], axis=plane).reshape(nshp)
		args = []
		if not coords is None:
			X = compress(cond, coords[(plane+1)%3], axis=plane).reshape(nshp)
			Y = compress(cond, coords[(plane-1)%3], axis=plane).reshape(nshp)
			args += [X[::dens[0],::dens[1]].T, Y[::dens[0],::dens[1]].T]
		args += [U[::dens[0],::dens[1]].T, V[::dens[0],::dens[1]].T]
		if not C is None:
			args += [C]
		quiver(*tuple(args), **kwargs)

def svprod(f, v):
	'''
Scalar-vector multiplication. The result is a vector (field).
	'''
	return VectorField([u*f for u in v])

def scalarprod(v, w):
	'''
Vector-vector scalar multiplication. The result is a scalar (field).
	'''
	if len(v) != len(w):
		raise TypeError('Invalid vectors')
	return sum([v[i]*w[i] for i in range(len(v))])

def coordCyl(*varargs):
	'''
Given three 1-D arrays representing cartesian the coordinates of a grid,
a centre and a direction, returns three vector fields containing the
cartesian coordinates of the basis vectors of a cylindric coordinate system.

Parameters
----------
x1, x2,..., xn : array_like
    1-D arrays representing the coordinates of a grid.
c : array_like (three components)
    centre of the new coordinate system
e : array_like (three components)
    vector in the direction of the axis of the (right-handed) cylinder
	'''
	if len(varargs)!=5:
		raise TypeError('You must provide three coordinates, a centre, and an axis')
	coords = varargs[:-2]
	c = varargs[-2]
	e = varargs[-1]
	if not len(c)==len(e)==3:
		raise TypeError('Invalid coordinate centre or cylinder axis')
	if e == 0.:
		raise TypeError('Cylinder axis cannot be zero')
	e0 = VectorField(asarray(e)/norm(e))
	c0 = VectorField(asarray(c))
	p = VectorField(mgrid(*coords))
	r  = (p-c0) - svprod(scalarprod(p-c0,e0), e0)
	theta = e0^r
	z = e0*ones_like(p)

	return [r, theta, z]
