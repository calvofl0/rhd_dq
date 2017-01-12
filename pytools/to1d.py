#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from pybold import level   # Slow routine coded in Python
from pylevel import level   # Fast routine coded in C
from numpy import mean, nanmean, asarray, arange, zeros_like, where, nan, isnan, empty, empty_like, shape, meshgrid, ceil
from scipy.interpolate.interpolate import spline
from numpy import count_nonzero, logical_not, argwhere, float32, asarray, zeros, newaxis, copyto
from numpy import min as npmin
from numpy import max as npmax
from _pybold import write_model

def zScale(tau_box, tau, isotau_flag=False):
    l = zeros_like(tau_box[:,:,0])
    z = []
    if isotau_flag:
        isotau = empty_like(tau_box)
    for i,t in zip(range(len(tau)),tau):
        print('Step: '+str(i))
        if isnan(t):
            z+=[nan]
            if isotau_flag:
                isotau[:,:,i] = nan
        else:
            l = level(tau_box, t, l)
            if isotau_flag:
                isotau[:,:,i] = l
            z += [nanmean(l)]

    if isotau_flag:
        return z, isotau
    return z

def tauScale(tau_box, isotau_flag=False):
    tau = mean(mean(tau_box,axis=0),axis=0)
    if isotau_flag:
        z, isotau = zScale(tau_box, tau, isotau_flag)
    else:
        z = zScale(tau_box, tau, isotau_flag)
    s = spline(z, tau, arange(len(z)))
    s[where(s==0.)] = nan

    if isotau_flag:
        return s, isotau
    return s

def meanVarAtIsotau(var, isotau):
    nx, ny, nz = shape(isotau)
    if type(var) is tuple:
        var1d = tuple([empty(nz) for i in range(len(var))])
    else:
        var1d = empty(nz)

    for i in range(nz):
        print('Step: '+str(i))
        nanvals = where(isnan(isotau[:,:,i]))
        frac = isotau[:,:,i]-ceil(isotau[:,:,i]-1.)
        isosurf = (isotau[:,:,i]-frac).astype(int)
        isosurf[nanvals] = 0
        #isosurf[where(isosurf>=nz-1)] = 0
        indices=tuple(list(meshgrid(range(ny),range(nx)))+[isosurf])
        indicesR=tuple(list(meshgrid(range(ny),range(nx)))+[isosurf+1])
        if type(var) is tuple:
            for j in range(len(var)):
                varAtIsotau = (1.-frac)*var[j][indices]+frac*var[j][indicesR]
                varAtIsotau[nanvals] = nan
                var1d[j][i] = nanmean(varAtIsotau)
        else:
            varAtIsotau = (1.-frac)*var[indices]+frac*var[indicesR]
            varAtIsotau[nanvals] = nan
            var1d[i] = nanmean(varAtIsotau)

    return var1d

def meanVarFlat(var):
    if type(var) is tuple:
        nz = shape(var[0])[2]
        var1d = tuple([empty(nz) for i in range(len(var))])
        for v,i in zip(var,range(nz)):
            copyto(var1d[i], mean(mean(v,axis=0),axis=0))
    else:
        var1d = mean(mean(var,axis=0),axis=0)

    return var1d


def export1d(model, outfile, nx=3, ny=3, flat=False, arrtype=float32):
    tau_box = model.dq.tau
    # Get ``equidistant'' tau-isosurfaces
    if not flat:
        tau, isotau = tauScale(tau_box, isotau_flag=True)
    # Average physical quantities on those surfaces (var is a tuple containing, for each quantity, a 1d-vector where each element is the average on the corresponding tau isosurface
    v1 = model.z['v1'].value
    v2 = model.z['v2'].value
    v3 = model.z['v3'].value
    rho = model.z['rho'].value
    ei = model.z['ei'].value
    if not flat:
        var = meanVarAtIsotau((v1, v2, v3, rho, ei), isotau)
    else:
        var = meanVarFlat((v1, v2, v3, rho, ei))
    # Cells where an isosurface was successfully computed
    if not flat:
        z = argwhere(logical_not(isnan(tau))).flatten()
    else:
        z = arange(rho.shape[2])
    m3 = z[0]
    n3 = z[-1]
    nz = n3-m3+1
    m1 = m2 = 1
    n1 = nx
    n2 = ny
    xc1 = model.z['xc1'].value.flatten()[arange(nx)]
    xc2 = model.z['xc2'].value.flatten()[arange(ny)]
    xc3 = model.z['xc3'].value.flatten()[z]
    xb1 = model.z['xb1'].value.flatten()[arange(nx+1)]
    xb2 = model.z['xb2'].value.flatten()[arange(ny+1)]
    xb3 = model.z['xb3'].value.flatten()[asarray(list(z)+[1+z[-1]])]
    v1 = asarray(var[0][z][newaxis,newaxis,:].repeat(nx,axis=0).repeat(ny,axis=1),order='F', dtype=arrtype)
    v2 = asarray(var[1][z][newaxis,newaxis,:].repeat(nx,axis=0).repeat(ny,axis=1),order='F', dtype=arrtype)
    v3 = asarray(var[2][z][newaxis,newaxis,:].repeat(nx,axis=0).repeat(ny,axis=1),order='F', dtype=arrtype)
    rho = asarray(var[3][z][newaxis,newaxis,:].repeat(nx,axis=0).repeat(ny,axis=1),order='F', dtype=arrtype)
    ei = asarray(var[4][z][newaxis,newaxis,:].repeat(nx,axis=0).repeat(ny,axis=1),order='F', dtype=arrtype)
    Bb_flag = False
    Bb1 = zeros((nx+1, ny, nz), dtype=arrtype, order='F')
    Bb2 = zeros((nx, ny+1, nz), dtype=arrtype, order='F')
    Bb3 = zeros((nx, ny, nz+1), dtype=arrtype, order='F')
    B1_unit = B2_unit = B3_unit = ''
    if hasattr(model.z,'time_db'):
            time_db=model.z['time_db'].value
    else: time_db=model['modeltime'].value
    if hasattr(model, 'head'):
            description=model.head['description'].value
            version=model.head['version'].value
            history=model.head['history'].value
    else :
            description=''
            version='Unknown version'
            history=np.array(['Unknown past history'])
    if hasattr(model, 'time_out_mean_last') :
            toml = model['time_out_mean_last'].value
    else :
            toml = 0.0
    if hasattr(model, 'time_out_full_last') :
            tofl = model['time_out_full_last'].value
    else :
            tofl = 0.0
    history=asarray([list(h.ljust(80)) for h in history])

    write_model(outfile, xb1, xb2, xb3, xc1, xc2, xc3, v1, v2, v3, rho, ei, Bb_flag, Bb1, Bb2, Bb3, B1_unit, B2_unit, B3_unit, model['dtime'].value, model['modelitime'].value, model['modeltime'].value, time_db, description, history, version, toml, tofl, m1, n1, m2, n2, m3, n3)
