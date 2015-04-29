#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
import numpy as np

def check_ext(filename) :
	if not filename.endswith('.out') :
		raise argparse.ArgumentTypeError("%s has not a recognized extension." % filename)
        return filename

description = 'Parse CO5BOLD out files to make statistics.'
parser = argparse.ArgumentParser(description=description)
parser.add_argument('input', help='CO5BOLD out file(s).', type=check_ext, nargs='+')
args = parser.parse_args()

phystime = []
cputime = []
mhdtime = []
rttime = []
visctime = []
uiotime = []
dtime = []
nitime = []
for filename in args.input :
	time0 = -1
	time1 = -1
	itime0 = -1
	itime1 = -1
	ldtime = []
	with open(filename) as f :
		lines = f.readlines()
	for l in lines :
		m = re.match(r'.*itime=    (.*)  time= (.*)  t_job',l)
		if m != None:
			if time0 < 0: time0 = float(m.group(2))
			if itime0 < 0: itime0 = float(m.group(1))
			time1 = float(m.group(2))
			itime1 = float(m.group(1))
		m = re.match(r'dtime= (.*) HD= (.*) RAD= (.*) VIS= (.*)',l)
		if m != None:
			ldtime += [float(m.group(1))]
		m = re.match(r'RHD code',l)
		if m != None:
			cputime += [float(l.split()[3])]
		m = re.match(r'Magneto-Hydrodynamics routines',l)
		if m != None:
			mhdtime += [float(l.split()[3])]
		m = re.match(r'Radiation transport routines',l)
		if m != None:
			rttime += [float(l.split()[4])]
		m = re.match(r'Viscosity routines 3D',l)
		if m != None:
			visctime += [float(l.split()[4])]
		m = re.match(r'uio output routines',l)
		if m != None:
			uiotime += [float(l.split()[4])]
	if time0 < 0: continue
	phystime += [time1-time0]
	nitime += [itime1-itime0]
	dtime += [np.mean(ldtime)]

phystime = np.array(phystime)
cputime = np.array(cputime)
mhdtime = np.array(mhdtime)
rttime = np.array(rttime)
dtime = np.array(dtime)
nitime = np.array(nitime)
rate = phystime/cputime*3600.*24.
mhdratio = mhdtime/cputime
rtratio = rttime/cputime
viscratio = visctime/cputime
uioratio = uiotime/cputime
tstep = cputime/nitime

tabs = 25
print('Simulation rate:\t{:0.0f} +/- {:0.0f} seconds/day'.format(np.mean(rate), np.std(rate)).expandtabs(tabs))
print('MHD time:\t{:0.0f}% +/- {:0.0f}%'.format(100.*np.mean(mhdratio), 100.*np.std(mhdratio)).expandtabs(tabs))
print('RT time:\t{:0.0f}% +/- {:0.0f}%'.format(100.*np.mean(rtratio), 100.*np.std(rtratio)).expandtabs(tabs))
print('Viscosity time:\t{:0.0f}% +/- {:0.0f}%'.format(100.*np.mean(viscratio), 100.*np.std(viscratio)).expandtabs(tabs))
print('I/O time:\t{:0.0f}% +/- {:0.0f}%'.format(100.*np.mean(uioratio), 100.*np.std(uioratio)).expandtabs(tabs))
print('CPU time per timestep:\t{:0.0f} +/- {:0.0f} seconds'.format(np.mean(tstep), np.std(tstep)).expandtabs(tabs))
