#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from array import array
from struct import pack as s_pack
from struct import unpack as s_unpack

def pack(tup, precision='d', output='chr'):
	if np.size(tup) == 1:
		try: len(tup)
		except: tup=(tup,)
	if precision not in set(['d', 'f']):
		raise ValueError("Optional argument precision must be either 'f' (single precision) or 'd' (double precision).")
	if output not in set(['chr', 'ord']):
		raise ValueError("Optional argument output must be wither 'chr' (pack into string) or 'ord' (pack into array of uint8)")
	l = []
	for var in tup:
		t   = None
		arr = False
		if type(var) == int:
			t='i'
		elif type(var) == float:
			t=precision
		elif type(var) == np.float32:
			t=precision
		elif type(var) == np.float64:
			t=precision
		elif type(var) == np.ndarray:
			arr = True
			if var.dtype == np.dtype(int):
				t='i'
			elif var.dtype == np.dtype(float):
				t=precision
			elif var.dtype == np.float32:
				t=precision
		if not t:
			if arr:
				raise ValueError("Unrecognized type: "+str(var.dtype)+".")
			else:
				raise ValueError("Unrecognized type: "+str(type(var))+".")
		if arr:
			l += [array(t, var).tostring()]
		else:
			l += [s_pack(t, var)]
	stream = ''.join(l)
	if output == 'chr':
		return stream
	else:
		out = np.empty((len(stream)), dtype=np.uint8)
		for i in range(len(stream)):
			out[i] = ord(stream[i])
		return out

def unpack(fmt, stream):
	if type(stream) == np.ndarray:
		stream0 = ''.join([chr(i) for i in stream])
	else:
		stream0 = stream
	return s_unpack(fmt, stream0)
