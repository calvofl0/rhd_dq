#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from binarytools import pack, unpack
from struct_type import Struct

NJKQ = 9

def readHeader(model):
	int_sz = model.int_sz
	db_sz  = model.db_sz

	offset = 0
	header = Struct()

	(header.version,) = unpack('i', model.module_header[offset:offset+int_sz])
	offset += int_sz
	if header.version != 1:
		raise Exception("Version of pratic module unsupported (only v1 supported).")
	(header.Nfreq,) = unpack('i', model.module_header[offset:offset+int_sz])
	offset += int_sz
	header.freq = np.array(list(unpack(str(header.Nfreq)+'d', model.module_header[offset:offset+db_sz*header.Nfreq])))
	offset += db_sz*header.Nfreq
	(header.nx, header.ny) = unpack('2i', model.module_header[offset:offset+db_sz])
	offset += db_sz
	header.x = np.array(list(unpack(str(header.nx)+'d', model.module_header[offset:offset+db_sz*header.nx])))
	offset += db_sz*header.nx
	header.y = np.array(list(unpack(str(header.ny)+'d', model.module_header[offset:offset+db_sz*header.ny])))
	offset += db_sz*header.ny
	header.T = np.array(list(unpack(str(header.ny*header.nx)+'d', model.module_header[offset:offset+db_sz*header.ny*header.nx]))).reshape((header.ny, header.nx)).T
	offset += db_sz*header.ny*header.nx

	return header

def writeHeader(header, model):
	h = []
	h += [pack(header.version)]
	header.Nfreq = len(header.freq)
	h += [pack(header.Nfreq)]
	h += [pack(header.freq)]
	(header.ny, header.nx) = (len(header.y), len(header.x))
	h += [pack(header.nx)]
	h += [pack(header.ny)]
	h += [pack(header.x)]
	h += [pack(header.y)]
	Tx, Ty = np.shape(header.T)
	if (Ty, Tx) != (header.ny, header.nx):
		raise Exception('Temperature has shape ('+Tx+','+Ty+') and expected shape is ('+header.nx+','+header.ny+').')
	h += [pack(header.T.T.flatten())]

	model.module_header = ''.join(h)

def readNodes(model):
	nx, ny, nz = model.dimensions
	node_data  = model.node_data
	header     = readHeader(model)
	Nfreq      = header.Nfreq
	int_sz     = model.int_sz
	db_sz      = model.db_sz

	node_data = Struct()
	node_data.ne      = np.empty(model.dimensions)
	node_data.nH1     = np.empty(model.dimensions)
	node_data.T       = np.empty(model.dimensions)
	node_data.kappa_c = np.empty((Nfreq, nx, ny, nz))
	node_data.JKQ     = np.empty((NJKQ, Nfreq, nx, ny, nz))
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				offset = 0
				(node_data.ne[i,j,k],) = unpack('d', model.node_data[i,j,k][offset:offset+db_sz])
				offset += db_sz
				(node_data.nH1[i,j,k],) = unpack('d', model.node_data[i,j,k][offset:offset+db_sz])
				offset += db_sz
				(node_data.T[i,j,k],) = unpack('d', model.node_data[i,j,k][offset:offset+db_sz])
				offset += db_sz
				for f in range(Nfreq):

					(node_data.kappa_c[f,i,j,k],) = unpack('d', model.node_data[i,j,k][offset:offset+db_sz])
					offset += db_sz
					node_data.JKQ[:,f,i,j,k] = unpack(str(NJKQ)+'d', model.node_data[i,j,k][offset:offset+db_sz*NJKQ])
					offset += db_sz*NJKQ
	return node_data

def writeNodes(node_data, model):
	nx, ny, nz = model.dimensions
	header     = readHeader(model)
	Nfreq      = header.Nfreq
	int_sz     = model.int_sz
	db_sz      = model.db_sz

	if not np.all(np.shape(node_data.ne) == model.dimensions):
		raise Exception("Variable 'ne' must have shape "+str(model.dimensions)+".")
	if not np.all(np.shape(node_data.nH1) == model.dimensions):
		raise Exception("Variable 'nH1' must have shape "+str(model.dimensions)+".")
	if not np.all(np.shape(node_data.T) == model.dimensions):
		raise Exception("Variable 'T' must have shape "+str(model.dimensions)+".")
	if not np.all(np.shape(node_data.kappa_c) == (Nfreq, nx, ny, nz)):
		raise Exception("Variable 'kappa_c' must have shape "+str((Nfreq,nx,ny,nz))+".")
	if not np.all(np.shape(node_data.JKQ) == (NJKQ, Nfreq, nx, ny, nz)):
		raise Exception("Variable 'JKQ' must have shape "+str((NJKQ,Nfreq,nx,ny,nz))+".")

	nd = np.empty((nx,ny,nz,db_sz*(3+(1+NJKQ)*Nfreq)), dtype=np.uint8)
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				offset = 0
				nd[i,j,k][:3*db_sz] = pack((node_data.ne[i,j,k], node_data.nH1[i,j,k], node_data.T[i,j,k]), output='ord')
				offset += 3*db_sz
				for f in range(Nfreq):
					nd[i,j,k][offset:offset+db_sz] = pack(node_data.kappa_c[f,i,j,k], output='ord')
					offset += db_sz
					nd[i,j,k][offset:offset+NJKQ*db_sz] = pack(node_data.JKQ[:,f,i,j,k], output='ord')
					offset += NJKQ*db_sz

	model.node_data = nd
