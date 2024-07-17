#!/usr/bin/env python

import optparse
import os,sys
from EMAN2 import *
import math


#===============
def setupParserOptions():
	parser = optparse.OptionParser()
	
	parser.set_usage("""%prog -f helical.mrc --DO --DI --DZ --boxsize --apix
			 to apply a soft mask to both ends of a 3D helical volume (z direction),
			 together with a soft inner and outer mask (x-y plane)""")
	parser.add_option("-f", dest="mrc", type="string", metavar="FILE",
		help="3D helical.mrc file to be masked. If not specified, will only output mask3D.mrc")
	parser.add_option("--DO", dest="DO", type="int", metavar="INT",
		help="outer diameter in x-y plane, in angstrom")
	parser.add_option("--DI", dest="DI", type="int", metavar="INT", default=0,
		help="inner diameter in x-y plane, in angstrom, default=0")
	parser.add_option("--DZ", dest="DZ", type="float", metavar="FLOAT",
		help="diameter in z direction, if >1, in angstrom; if <1, in percentage")
	parser.add_option("--boxsize", dest="boxsize", type="int", metavar="INT",
		help="boxsize")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT",
		help="apix")
	
	options,args = parser.parse_args()
	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))
	if len(sys.argv) < 2:
		parser.print_help()
		parser.error("no options defined")
	params = {}
	for i in parser.option_list:
		if isinstance(i.dest, str):
			params[i.dest] = getattr(options, i.dest)
	return params

# Generate 2D slices to be inserted into mask3D volume
def createMask2D(params):
	from itertools import product
	apix = params['apix']
	RO = float(params['DO'])/apix/2
	RI = float(params['DI'])/apix/2
	if params['boxsize']:
		nx = int(params['boxsize'])
	else:
		vol = EMData(params['mrc'])
		nx = vol.get_xsize()
	falloff_r = 6			# use steeper falloff
	mask2D = EMData(nx,nx)
	mask2D.to_one()
	for x,y in product(range(nx),range(nx)):
		dx = abs(x-nx/2)
		dy = abs(y-nx/2)
		r2 = dx**2+dy**2
		if r2 > RO*RO:
			wt1 = 0.5*(1 + math.cos(math.pi*min(1,(math.sqrt(r2)-RO)/falloff_r)))
			mask2D.set(x,y,wt1)
		elif r2 < RI*RI:
			wt2 = 0.5*(1 + math.cos(math.pi*min(1,(RI-math.sqrt(r2))/falloff_r)))
			mask2D.set(x,y,wt2)
	#mask2D.write_image('mask.mrc')
	return mask2D

def createMask3D(params,mask2D):
	apix = params['apix']
	RO = float(params['DO'])/apix/2
	RI = float(params['DI'])/apix/2
	nx = mask2D.get_xsize()
	if params['DZ'] > 1.0:
		RZ = float(params['DZ'])/apix/2
	else:
		RZ = (nx/2*params['DZ'])
	mask3D = EMData(nx,nx,nx)
	falloff_z = 6.0
	# now apply soft mask
	for z in range(nx):
		img = EMData(nx,nx)
		img = mask2D.copy()
		# here "img = mask2D" won't work !!
		dz = abs(z-nx/2)
		if dz > RZ:
			wt3 = 0.5*(1 + math.cos(math.pi*min(1,(dz-RZ)/falloff_z)))
			img.mult(wt3)
			#img.write_image("test_%d.mrc"%z)
		mask3D.insert_clip(img,(0,0,z))
	
	name = 'mask3D_%d-%dx%d_apix%.2f.mrc'%(params['DO'],params['DI'],math.ceil(RZ*apix*2),apix)
	mask3D.write_image(name)
	return mask3D,name




#===============
if  __name__ == "__main__":
	params = setupParserOptions()
	mask2D = createMask2D(params)
	mask3D,name = createMask3D(params,mask2D)	
	if params['mrc']:
		vol = EMData(params['mrc'])
		vol.mult(mask3D)
		vol.write_image('%s_mask.mrc'%(params["mrc"][:-4]))
	#os.system("bhead -origin 0,0,0 -recalculate -sampling %.3f %s %s"%(params['apix'],name,name))
