#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python interface to the pyvisit backend Fortran module.

See documentation for the "simulate" function.
"""

from _pyvisit import *
import pybold
from numpy import shape
from numpy import size
from numpy import ndarray
import numpy as np

arrtype = np.float32

# Testing

#mean=pybold.uio_struct()
#mean.load('rhd-test.mean')
#mean.objects
#model=pybold.uio_struct()
#model.load('rhd-test.sta')
#model.objects

#pv_regsim("TestSim", "Test simulation", "/no/useful/path")
#pv_regmesh2d("mesh2d", mean.rad['xb1'].value, mean.rad['xb2'].value, mean.rad['xb1'].unit, mean.rad['xb2'].unit)
#pv_regscalar("rad.intb3_r", "mesh2d", mean.rad['intb3_r'].value, mean.rad['intb3_r'].unit)
#pv_runsim()

def simulate(*args, **kwargs):
	'''
	Lauch a simulation with VisIt.

	>>> pyvisit.simulate(var1, (field1_x, field1_y, field1_z), var2,...)

	One can also provide optional arguments:
		- max_sz	(maximal size of the box sent to VisIt, make it
				coarser if necessary)
		- box   	(the box is first of all cropped to this size)
		- centre 	(when cropping, specify the centre of the new
				smaller box).
	'''
	meshes = ()	# Tuple of uio_struct (all different) corresponding to args
	mesh_id = ()	# Tuple of indices, so that meshes[mesh_id[i]] corresponds to args[i]
	dim = ()	# Tuple of meshes dimensions for each arg in args
	mesh_name = ()	# ............... names      ....................
	mesh_nodes = () # ............... #nodes     ....................
	mesh_cells = () # ............... #cells     ....................
	mesh_step = ()	# ............... steps      ....................
			# The new coarser mesh will only have one cell out
			# of steps. Each step is a 3-d tuple for each dimension
			# No coarser grid in 2d
	trash = ()	# Contains variables that are no more used in Python but
			# that need to stay in memory for Fortran routines
	# Parse optional arguments
	for key in kwargs:
		if key=='max_sz': max_sz=kwargs[key]
		elif key=='box': box=kwargs[key]
		elif key=='centre': centre=kwargs[key]
	if not 'max_sz' in locals(): max_sz=(-1,-1,-1)
	if not 'box' in locals(): box=(-1,-1,-1)
	if not 'centre' in locals(): centre=(-1,-1,-1)
	def tuple_swap_xy(t):
		u=(t[1], t[0])
		u+=t[2:]
		return u
	max_sz=tuple_swap_xy(max_sz)
	box=tuple_swap_xy(box)
	centre=tuple_swap_xy(centre)
	
	# Parse boxes to plot and identify relevant meshes
	index = 0
	for key in args:
		if type(key) == tuple:
			key=key[0]
		if type(key) == pybold.uio_struct_item:
			key = key._parent
		has_xb=False
		has_dq=False
		has_xb=hasattr(key,'xb1') and hasattr(key,'xb2') and hasattr(key,'xb3')
		if(hasattr(key,'box_id')):
			if(key['box_id'].value=='dq'):
				has_dq=hasattr(key._parent.z,'xb1') and hasattr(key._parent.z,'xb2') and hasattr(key._parent.z,'xb3')
		if(not key in meshes):
			if(has_xb):
				meshes += (key,)
			elif(has_dq):
				meshes += (key._parent.z,)
		if(has_xb):
			mesh_id += (np.where(np.array(meshes)==key)[0][0].item(),)
		elif(has_dq):
			mesh_id += (np.where(np.array(meshes)==key._parent.z)[0][0].item(),)
		else:
			mesh_id += (None,)
		index += 1
	# Initialize new simulation
	pv_reset()
	pv_regsim("rhd_dq_sim", "CO5BOLD box simulation", "/no/useful/path")
	# Register new meshes
	index=0
	for mesh in meshes:
		if(mesh['box_id'].value=='rad'): dim+=(2,)
		elif(shape(mesh['xb1'].value)[0]<=2 or shape(mesh['xb2'].value)[1]<=2 or shape(mesh['xb3'].value)[2]<=2): dim+=(2,)
		else: dim+=(3,)
		mesh_name += (mesh['box_id'].value+'.mesh'+str(mesh['itime'].value),)
		if(dim[index]==2):
			if shape(mesh['xb3'].value)[2]<=2:
				pv_regmesh2d(mesh_name[index],mesh['xb1'].value,mesh['xb2'].value,mesh['xb1'].unit,mesh['xb2'].unit)
				mesh_nodes += (shape(mesh['xb1'].value)[0]*shape(mesh['xb2'].value)[1],)
				mesh_cells += ((shape(mesh['xb1'].value)[0]-1)*(shape(mesh['xb2'].value)[1]-1),)
			elif shape(mesh['xb2'].value)[1]<=2:
				pv_regmesh2d(mesh_name[index],mesh['xb1'].value,mesh['xb3'].value,mesh['xb1'].unit,mesh['xb3'].unit)
				mesh_nodes += (shape(mesh['xb1'].value)[0]*shape(mesh['xb3'].value)[2],)
				mesh_cells += ((shape(mesh['xb1'].value)[0]-1)*(shape(mesh['xb3'].value)[2]-1),)
			elif shape(mesh['xb1'].value)[0]<=2:
				pv_regmesh2d(mesh_name[index],mesh['xb2'].value,mesh['xb3'].value,mesh['xb2'].unit,mesh['xb3'].unit)
				mesh_nodes += (shape(mesh['xb2'].value)[1]*shape(mesh['xb3'].value)[2],)
				mesh_cells += ((shape(mesh['xb2'].value)[1]-1)*(shape(mesh['xb3'].value)[2]-1),)
			else:
				pv_regmesh2d(mesh_name[index],mesh['xb1'].value,mesh['xb2'].value,mesh['xb1'].unit,mesh['xb2'].unit)
				mesh_nodes += (shape(mesh['xb1'].value)[0]*shape(mesh['xb2'].value)[1],)
				mesh_cells += ((shape(mesh['xb1'].value)[0]-1)*(shape(mesh['xb2'].value)[1]-1),)
			mesh_step += (None,)
		else:
			# 3D meshes need to be made coarser
			id_min0=np.array([0,0,0])
			id_max0=np.array([shape(mesh['xb1'].value)[0],shape(mesh['xb2'].value)[1],shape(mesh['xb3'].value)[2]])
			id_min=id_min0.copy()
			id_max=id_max0.copy()
			steps=np.array([1,1,1])
			for i in range(3):
				if centre[i]<0 and box[i]>=0 : centre[i]=0
				if centre[i]>=0 :
					if box[i]>=0 :
						id_min[i]=max(id_min0[i],int(np.ceil(centre[i]-box[i]/2.)))
						id_max[i]=min(id_max0[i],1+int(np.floor(centre[i]+box[i]/2.)))
						if max_sz[i]>=0 :
							steps[i]=max(1,int(np.floor((id_max[i]-id_min[i]+1)/max_sz[i])))
					else :
						if max_sz[i]>=0 :
							id_min[i]=max(id_min[i],int(np.ceil(centre[i]-max_sz[i]/2.)))
							id_max[i]=min(id_max[i],1+int(np.floor(centre[i]+max_sz[i]/2.)))
				else :
					if max_sz[i]>=0 :
						steps[i]=max(1,int(np.floor((id_max[i]-id_min[i]+1)/max_sz[i])))
			if id_min0[0]==id_min[0] and id_max0[0]==id_max[0] and steps[i]==1 :
				xb1=np.asfortranarray(mesh['xb1'].value)
			else :
				xb1=np.array(mesh['xb1'].value[id_min[0]:id_max[0]:steps[0],:,:],order='F', dtype=arrtype)
			trash+=(xb1,)
			if id_min0[1]==id_min[1] and id_max0[1]==id_max[1] and steps[i]==1 :
				xb2=np.asfortranarray(mesh['xb2'].value)
			else :
				xb2=np.array(mesh['xb2'].value[:,id_min[1]:id_max[1]:steps[1],:],order='F', dtype=arrtype)
			trash+=(xb2,)
			if id_min0[2]==id_min[2] and id_max0[2]==id_max[2] and steps[i]==1 :
				xb3=np.asfortranarray(mesh['xb3'].value)
			else :
				xb3=np.array(mesh['xb3'].value[:,:,id_min[2]:id_max[2]:steps[2]],order='F', dtype=arrtype)
			trash+=(xb3,)
			pv_regmesh3d(mesh_name[index],xb2,xb1,xb3,mesh['xb2'].unit,mesh['xb1'].unit,mesh['xb3'].unit)
			mesh_nodes += (shape(mesh['xb1'].value)[0]*shape(mesh['xb2'].value)[1]*shape(mesh['xb3'].value)[2],)
			mesh_cells += ((shape(mesh['xb1'].value)[0]-1)*(shape(mesh['xb2'].value)[1]-1)*(shape(mesh['xb3'].value)[2]-1),)
			mesh_step += (steps,)
		index += 1
	# Expand arguments in args, so that an uio_struct argument will be
	# expanded in a list of uio_struct_item arguments and a 3d-tuple of
	# uio_struc_item will stay unchanged
	index = 0
	values = ()
	mesh_id0 = ()
	for key in args:
		if type(key) == pybold.uio_struct:
			for item in key.viewkeys():
				if type(key[item]) == pybold.uio_struct_item:
					if type(key[item].value) == ndarray and not item in values:
						values += (key[item],)
						mesh_id0 += (mesh_id[index],)
		elif type(key) == pybold.uio_struct_item:
			values += (key,)
			mesh_id0 += (mesh_id[index],)
		elif type(key) == tuple:
			if len(key) == 3:
				if type(key[0])==pybold.uio_struct_item and type(key[1])==pybold.uio_struct_item and type(key[2])==pybold.uio_struct_item:
					values += (key,)
					mesh_id0 += (mesh_id[index],)
		index += 1
	# Register all scalars and vectors
	index = 0
	for key in values:
		if type(key) == tuple:
			if np.all(np.prod([np.array(np.shape(key[0].value))-np.array([1,0,0]),np.array(np.shape(key[1].value))-np.array([0,1,0]),np.array(np.shape(key[2].value))-np.array([0,0,1])],axis=1) == mesh_cells[mesh_id0[index]]):
				box1 = .5*(key[0].value[0:-1,:,:]+key[0].value[1:,:,:])
				box2 = .5*(key[1].value[:,0:-1,:]+key[1].value[:,1:,:])
				box3 = .5*(key[2].value[:,:,0:-1]+key[2].value[:,:,1:])
				print('Reg b-vect')
			elif np.all(np.array([size(key[0].value), size(key[1].value), size(key[2].value)]) == mesh_cells[mesh_id0[index]]):
				box1 = key[0].value
				box2 = key[1].value
				box3 = key[2].value
				print('Reg vect')
			else: continue
			#print(str(index))
			#print(str(mesh_id0[index]))
			#print(mesh_step)
			#print(mesh_step[mesh_id0[index]])
			box1 = box1[id_min[0]:id_max[0]-1:mesh_step[mesh_id0[index]][0],id_min[1]:id_max[1]-1:mesh_step[mesh_id0[index]][1],id_min[2]:id_max[2]-1:mesh_step[mesh_id0[index]][2]]
			box2 = box2[id_min[0]:id_max[0]-1:mesh_step[mesh_id0[index]][0],id_min[1]:id_max[1]-1:mesh_step[mesh_id0[index]][1],id_min[2]:id_max[2]-1:mesh_step[mesh_id0[index]][2]]
			box3 = box3[id_min[0]:id_max[0]-1:mesh_step[mesh_id0[index]][0],id_min[1]:id_max[1]-1:mesh_step[mesh_id0[index]][1],id_min[2]:id_max[2]-1:mesh_step[mesh_id0[index]][2]]
			box=np.empty(np.size(box1)+np.size(box2)+np.size(box3), order='F', dtype=arrtype)
			box[0::3] = np.transpose(box1,(1,0,2)).flatten()
			box[1::3] = np.transpose(box2,(1,0,2)).flatten()
			box[2::3] = np.transpose(box3,(1,0,2)).flatten()
			trash+=(box,)
			pv_regvector(meshes[mesh_id0[index]]['box_id'].value+'.vect,'+key[0].name+','+key[1].name+','+key[2].name,mesh_name[mesh_id0[index]],box,'vect,'+key[0].unit+','+key[1].unit+','+key[2].unit)
		elif size(key.value)==mesh_cells[mesh_id0[index]]:
			if dim[mesh_id0[index]] == 2:
				box=np.array(key.value.T, order='F', dtype=arrtype)
				trash+=(box,)
			elif dim[mesh_id0[index]] == 3:
				box=np.array(np.transpose(key.value[id_min[0]:id_max[0]-1:mesh_step[mesh_id0[index]][0],id_min[1]:id_max[1]-1:mesh_step[mesh_id0[index]][1],id_min[2]:id_max[2]-1:mesh_step[mesh_id0[index]][2]],(1,0,2)), order='F', dtype=arrtype)
				trash+=(box,)
			print(meshes[mesh_id0[index]]['box_id'].value+'.'+key.name)
			print(mesh_name[mesh_id0[index]])
			pv_regscalar(meshes[mesh_id0[index]]['box_id'].value+'.'+key.name,mesh_name[mesh_id0[index]],box,key.unit)
		index += 1
	trash+=(mesh_cells, mesh_id)
	return trash

def simulate2(model, values=""):
	if values == "":
		values0=()
		for k in model.viewkeys():
			if type(model[k]==pybold.uio_struct):
				values0+=(k,)
	else:
		values0=tuple(values.split(','))
	items=()
	for val in values0:
		if type(model[val])==pybold.uio_struct:
			for k in model[val].viewkeys():
				items+=(val+'.'+k,)
	boxes={}
	for q in items:
		q=q.split('.')
		if not boxes.has_key(q[0]):
			boxes[q[0]]=[]
		boxes[q[0]]+=[q[1]]

	return boxes
