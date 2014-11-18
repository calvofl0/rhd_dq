#!/usr/bin/env python
# -*- coding: utf-8 -*-

PYBOLD=True
PIL=True
try:
	import pybold
except:
	PYBOLD=False
try:
	from PIL import Image
except:
	PIL=False
import argparse
import sys
import pickle
from numpy import min as npmin
from numpy import max as npmax
from matplotlib.pyplot import imsave
from matplotlib.pylab import cm
from tempfile import mkdtemp
from shutil import rmtree
from os import system
from os.path import basename
from os.path import abspath
from glob import glob
from socket import gethostname

digits = 3

def check_ext_in(filename) :
	if not (filename.endswith('.tar') or filename.endswith('.mean')) :
		raise argparse.ArgumentTypeError("%s has not a recognized extension." % filename)
	if filename.endswith('.mean') and not PYBOLD :
		raise argparse.ArgumentTypeError("You need pybold to open %s." % filename)
	return filename

def check_ext_out(filename) :
	if not (filename.endswith('.tar') or filename.endswith('.mp4')) :
		raise argparse.ArgumentTypeError("%s has not a recognized extension." % filename)
	return filename

def check_legend(legend) :
	if legend and not PIL :
		raise argparse.ArgumentTypeError("You need PIL Python module to add a legend.")
	return legend

description = 'Creates a movie from a MEAN file.'
parser = argparse.ArgumentParser(description=description)
parser.add_argument('-l', '--legend', help='add legend', action='store_true')
parser.add_argument('-m', '--marks', help='draw marks')
parser.add_argument('input', help='input MEAN/tar file.', type=check_ext_in, nargs='+')
parser.add_argument('output', help='output mp4/tar file.', type=check_ext_out)
args = parser.parse_args()

check_legend(args.legend)

tmpfolder = mkdtemp()

if args.input[0].endswith('.tar') :
	print('Extracting files...')
	system('tar xf '+args.input[0]+' -C '+tmpfolder)
	with open(tmpfolder+'/minmax.dat', 'rb') as f :
		vmin, vmax = pickle.load(f)
	if args.legend or args.marks :
		print('Adding legend and/or marks...')
		if args.marks :
			mrk_cdef  = 'red'
			mrk_rad   = 25
			mrk_space = 5
			mrk_label = []
			mrk_posx  = []
			mrk_posy  = []
			mrk_snap0 = []
			mrk_snap1 = []
			mrk_color = []
			mrk_serie = []
			with open(args.marks, 'rt') as f :
				marks=f.read().splitlines()
				serie=0
				for mark in marks :
					if mark.strip()=='': continue
					if mark.lstrip().startswith('#'): continue
					if mark.strip().startswith('---') :
						serie += 1
						continue
					mark=mark.split('\t')
					mrk_label += [mark[0]]
					mrk_posx  += [mark[1]]
					mrk_posy  += [mark[2]]
					mrk_snap0 += [mark[3]]
					mrk_snap1 += [mark[4]]
					mrk_serie += [serie]
					if len(mark) > 5 :
						mrk_color += [mark[5]]
					else : mrk_color += [mrk_cdef]
		with open(tmpfolder+'/index.dat', 'rt') as f :
			findex=f.read().splitlines()
		snpidx=0
		serie=0
		itime_prev=0
		for line in findex[1:] :
			entry=line.split('\t')
			print('\t'+entry[3])
			snap=Image.open(tmpfolder+'/'+entry[3])
			width,height=snap.size
			cmd_snap = entry[3]
			cmd_sidx = entry[2]
			cmd_mean = entry[1]
			cmd_time = str(round(10.*float(entry[5]))/10.)
			cmd_itime = entry[4]
			cmd_file = tmpfolder+'/'+cmd_snap
			cmd_leg  = ''
			if itime_prev > int(cmd_itime) : serie += 1
			itime_prev = int(cmd_itime)
			if args.legend :
				cmd_leg += ' \( \( -gravity NorthWest'
				cmd_leg += ' -background "#0008"'
				cmd_leg += ' -pointsize 24'
				cmd_leg += ' -size '+str(width)+'x30'
				cmd_leg += ' -fill "#ffff"'
				cmd_leg += ' caption:"'+cmd_snap+'"'
				cmd_leg += ' -background "#0000"'
				cmd_leg += ' -gravity NorthEast'
				cmd_leg += ' caption:"t='+cmd_time+'('+cmd_itime+')"'
				cmd_leg += ' -compose Over'
				cmd_leg += ' -composite'
				cmd_leg += ' -alpha off \)'
				cmd_leg += ' \( -gravity NorthWest'
				cmd_leg += ' -background "#0008"'
				cmd_leg += ' -pointsize 14'
				cmd_leg += ' -size '+str(width)+'x30'
				cmd_leg += ' -fill "#ffff"'
				cmd_leg += ' caption:"'+cmd_mean+'"'
				cmd_leg += ' -background "#0000"'
				cmd_leg += ' -gravity NorthEast'
				cmd_leg += ' caption:s="'+cmd_sidx+'"'
				cmd_leg += ' -compose Over'
				cmd_leg += ' -composite'
				cmd_leg += ' -alpha off \)'
				cmd_leg += ' -append \)'
				cmd_leg += ' -append'
			cmd_mrk = ''
			if args.marks :
				for i in range(len(mrk_label)) :
					if int(mrk_snap0[i]) <= int(cmd_itime) <= int(mrk_snap1[i]) and mrk_serie[i]==serie :
						cmd_mrk += ' -gravity NorthWest'
						cmd_mrk += ' -stroke '+mrk_color[i]
						cmd_mrk += ' -fill none'
						cmd_mrk += ' -draw "circle '+mrk_posx[i]+','+mrk_posy[i]+' '+str(int(mrk_posx[i])+mrk_rad)+','+mrk_posy[i]+'"'
						cmd_mrk += ' -draw "circle '+str(int(mrk_posx[i])+width)+','+mrk_posy[i]+' '+str(int(mrk_posx[i])+width+mrk_rad)+','+mrk_posy[i]+'"'
						cmd_mrk += ' -draw "circle '+mrk_posx[i]+','+str(int(mrk_posy[i])+height)+' '+str(int(mrk_posx[i])+mrk_rad)+','+str(int(mrk_posy[i])+height)+'"'
						cmd_mrk += ' -draw "circle '+str(int(mrk_posx[i])+width)+','+str(int(mrk_posy[i])+height)+' '+str(int(mrk_posx[i])+width+mrk_rad)+','+str(int(mrk_posy[i])+height)+'"'
						cmd_mrk += ' -stroke none'
						cmd_mrk += ' -fill '+mrk_color[i]
						cmd_mrk += ' -pointsize 20'
						cmd_mrk += ' -draw "text '+str(int(mrk_posx[i])+mrk_rad+mrk_space)+','+mrk_posy[i]+' \''+mrk_label[i]+'\''+'"'
						cmd_mrk += ' -draw "text '+str(int(mrk_posx[i])+width+mrk_rad+mrk_space)+','+mrk_posy[i]+' \''+mrk_label[i]+'\''+'"'
						cmd_mrk += ' -draw "text '+str(int(mrk_posx[i])+mrk_rad+mrk_space)+','+str(int(mrk_posy[i])+height)+' \''+mrk_label[i]+'\''+'"'
						cmd_mrk += ' -draw "text '+str(int(mrk_posx[i])+width+mrk_rad+mrk_space)+','+str(int(mrk_posy[i])+height)+' \''+mrk_label[i]+'\''+'"'
			cmd_action = ' -define png:color-type=2'
			system('convert '+cmd_file+cmd_mrk+cmd_leg+cmd_action+' '+cmd_file)
			snpidx+=1
			#system('convert '+tmpfolder+'/'+entry[3]+' \( -gravity NorthWest -background "#0008" -pointsize 24 -size '+str(width)+'x72 -fill "#ffff" caption:"'+entry[3]+'\n'+entry[1]+'" -background "#0000" -gravity NorthEast caption:"t='+entry[5]+'('+entry[4]+')'+'" -compose Over -composite -alpha off \) -append -define png:color-type=2 '+tmpfolder+'/'+entry[3])
elif args.input[0].endswith('.mean') :
	print("Normalizing boxes...")
	findex = open(tmpfolder+'/index.dat', 'wt')
	findex.write('Hostname\tMEAN file\tSnapshot index\tPNG file\titime\ttime\n')
	data = pybold.uio_struct()
	vmin = []
	vmax = []
	itime = []
	time = []
	selected = []
	nmodel = []
	itime_prev = 0
	for j in range(len(args.input)):
		i = 0
		data.open(args.input[j])
		while data.next > 0 :
			if data.modelitime > itime_prev or i==0 :
				vmin += [npmin(data.rad.intb3_r[:,:,0])]
				vmax += [npmax(data.rad.intb3_r[:,:,0])]
				selected += [True]
				itime_prev=int(data.modelitime)
				print('\t'+basename(args.input[j])+'\t'+str(data.modeltime)+'('+str(data.modelitime)+')')
			else :
				selected += [False]
			itime += [int(data.modelitime)]
			time += [float(data.modeltime)]
			i += 1
		data.close
		nmodel += [i]

	vmin = min(vmin)
	vmax = max(vmax)
	mean_index = 0
	print('Processing MEAN files...')
	all_index=0
	for j in range(len(args.input)):
		data.open(args.input[j])
		for i in range(nmodel[j]) :
			data.next
			if selected[all_index] :
				print('\t'+basename(args.input[j])+'\t'+str(time[all_index])+'('+str(itime[all_index])+')\t-->\t'+'snapshot'+('{0:0>'+str(digits)+'}').format(mean_index)+'.png')
				imsave(tmpfolder+'/snapshot'+('{0:0>'+str(digits)+'}').format(mean_index)+'.png',data.rad.intb3_r[:,:,0],vmin=vmin,vmax=vmax,cmap=cm.gray, origin='lower')
				findex.write(gethostname()+'\t'+abspath(args.input[j])+'\t'+('{0:0>'+str(digits)+'}').format(i+1)+'\t'+'snapshot'+('{0:0>'+str(digits)+'}').format(mean_index)+'.png\t'+str(itime[all_index])+'\t'+str(time[all_index])+'\n')
				mean_index += 1
			all_index+=1
		data.close
	findex.close()

if args.output.endswith('.tar') :
	with open(tmpfolder+'/minmax.dat', 'wb') as f :
		pickle.dump([vmin, vmax], f, protocol=-1)
	print('Creating TAR file...')
	system('tar cf '+args.output+' -C '+tmpfolder+' '+" ".join(map(lambda file:'"'+basename(file)+'"',glob(tmpfolder+'/*'))))
elif args.output.endswith('.mp4') :
	print('Creating movie...')
	system('ffmpeg -r 15 -i '+tmpfolder+'/snapshot%03d.png -c:v libx264 '+args.output)

rmtree(tmpfolder)
