#!/usr/bin/env python

from EMAN2 import *
import math
import matplotlib.pyplot as plt
import optparse
import os,sys
import numpy as np

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -f <EB3clK_5.par>")
        #parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
        #        help="EB3clK_5.par")
	parser.add_option("-s",dest="stack",type="string",metavar="FILE",
		help="start.hed")
	parser.add_option("-v",dest="vol",type="string",metavar="FILE",
		help="refvol.mrc")
	parser.add_option("--oversamp",dest="oversamp",type="int",default=2,
		help="oversampling factor")
	parser.add_option("--apix",dest="apix",type="float",default=1.32,
		help="apix")
	parser.add_option("--res",dest="res",type="float", default=4.8,
                help="use data up to certain resolution, default = 4.8")
	parser.add_option("--savedata", action="store_true",dest="savedata",default=True,
                help="savedata")
	
        options,args = parser.parse_args()

        if len(args) > 1:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit()

        params={}

        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0:
       return v
    return v/norm

def get1Dprofile(img,oversamp,Rad):
	# img need to be rotated to horizontal position
	# pad the image 2X
	nx = img.get_xsize()
	ny = img.get_ysize()
	nyo = ny*oversamp

	img.clip_inplace(Region(-(ny*(oversamp-1)/2),-(ny*(oversamp-1)/2),nyo,nyo))

	f = img.do_fft()
	f.process_inplace("xform.fourierorigin.tocenter")
	f.ri2ap()
	amp = f.amplitude()
	#display(amp)

	# calculate 1D profile
	nx = amp.get_xsize()
	ny = amp.get_ysize()

	# due to oversamp = 2
	Rad *= 2

	values = []
	for i in range(0,Rad):
		values.append(0)

	
	for x in range(0,Rad):
		Ry = int(math.sqrt(Rad**2-x**2))
		count = 0
		for y in range(ny/2-Ry,ny/2+Ry+1):
			values[x] += amp.get_value_at(x,y)
			count += 1
		#print x,count
		values[x] = values[x]/count
	return np.array(values)
	#return normalize(np.array(values))

def find2peaks(loc,mylist):
	# 7 pixel window
	mylist2 = list(mylist)
	mysublist = []
	for i in range(-3,4):
		mysublist.append(mylist2[loc+i])
	peak_1st = max(mysublist)
	index_1st = mylist2.index(peak_1st)
	peak_2nd = max(mylist2[index_1st-1],mylist2[index_1st+1])
	index_2nd = mylist2.index(peak_2nd)
	return peak_1st,peak_2nd,index_1st,index_2nd

def mainloop(params):
	oversamp = params['oversamp']
	apix = params['apix']
	res = params['res']
	if params['vol']:
		vol = EMData(params['vol'])
		Tp = Transform({"type":"mrc","theta":90.0,"phi":0.0})
		img = vol.project("standard",Tp)
		name = params['vol'][:-4]
	elif params['stack']:
		img = EMData(params['stack'],0)
		name = params['stack'][:-4]
	else:
		print "please specify input"
		sys.exit()
	if params['savedata']:
		if params['vol']:
			fout3 = file('%s.peak'%params['vol'][:-4],"w")
		elif params['stack']:
			fout3 = file('%s.peak'%params['stack'][:-4],"w")
	#display(img)
	nx = img.get_xsize()
	#ny = img.get_ysize()
	#a.process_inplace("normalize.edgemean")
	t1 = Transform()
	#psi = 90.0
	psi = 0.0
	t1.set_rotation({"type":"2d","alpha":-psi})
	img.process_inplace("xform",{"transform":t1})
	
	Rad = int(nx*apix/res)

	values = get1Dprofile(img,oversamp,Rad)
	# 512*1.32*2/4 = 270 pixels
	ub = int(512*1.32*2/res)
	background = min(values[2:ub])
	values -= background
	if params['vol']:
		fout2 = file('%s.1dp'%params['vol'][:-4],"w")
	elif params['stack']:
		fout2 = file('%s.1dp'%params['stack'][:-4],"w")

	if params['savedata']:
		# has to pre-calculate the peak40A_1st to scale up
		loc_40A = int(nx*apix/41*oversamp)
		loc_80A = int(nx*apix/82*oversamp)
		#background = min(values[loc_80A-2:loc_40A+2])
		#values -= background
		peak40A_1st,peak40A_2nd,index40A_1st,index40A_2nd = find2peaks(loc_40A,values)
		peak80A_1st,peak80A_2nd,index80A_1st,index80A_2nd = find2peaks(loc_80A,values)
		# sum the highest 2 values from the 5-values-window
		peaks80A = peak80A_1st + peak80A_2nd
		peaks40A = peak40A_1st + peak40A_2nd
		ratio80to40 = peak80A_1st/peak40A_1st
		ratios80to40 = peaks80A/peaks40A
		fout3.write("%.1f  %.1f  %.1f  %.1f\n"%(peak80A_1st,peak80A_2nd,peak40A_1st,peak40A_2nd))
		fout3.write("1 peak ratio 80/40 = %.2f\n"%ratio80to40)
		fout3.write("2 peaks ratio 80/40 = %.2f\n"%ratios80to40)
		fout3.close()

	values = values/peak40A_1st
	for x in range(0,Rad*2-1):
                #logvals.append(math.log(values[x]))
                fout2.write("%.6f\n"%values[x])
        fout2.close()

	#plt.plot(range(len(values)-2),values[2:])
	#plt.show()


if __name__ == "__main__":
        params = setupParserOptions()
	mainloop(params)

