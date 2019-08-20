#!/usr/bin/env python

import os,sys
import optparse
from EMAN2 import *
import math

#===============
def setupParserOptions():
	parser = optparse.OptionParser()
	
	#parser.set_usage("%prog -f fscData1.txt,fscData2.txt -m eman2/frealign --apix --boxsize --showres")
	parser.add_option("-f", dest="fnames", type="string",
		help="fsc data text files")
	parser.add_option("--boxsize", dest="boxsize", type="int", metavar="INT", default=512,
		help="box size, default=512")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT", default=1.403,
		help="pixel size in angstroms")
	parser.add_option("--res",dest="res",type="float", default=5.0,
                help="use data up to certain resolution, default = 5.0")
	parser.add_option("--drawpeak", action="store_true",dest="drawpeak", default=False,
                help="draw a line at the peak")
	
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

def parseDataFiles(params):
	fnames = params['fnames']
	fnameslist = fnames.split(',')
	
	#apix = params['apix']
	#boxsize = params['boxsize']
	
	# resList = [[res1],[res2],[res3]]
	ampList = []
	
	for i in fnameslist:
		f1 = file(i)
		l1 = f1.readlines()
		#amp = [float(i.strip())/1e6 for i in l1]
		#amp = [float(i.strip())/1e4 for i in l1]
		amp = [float(i.strip()) for i in l1]
		ampList.append(amp)
		f1.close()
	return ampList

def plotData(ampList,params):
	import matplotlib
	import matplotlib.pyplot as plt
	
	oversamp = 2
	# here apix = apix*2 due to oversampling
	apix = float(params['apix'])*oversamp
	boxsize = params['boxsize']
	print "boxsize = %d"%boxsize
	print "apix = %.2f"%apix

	loc_40A = int(boxsize*apix/41)
	loc_80A = int(boxsize*apix/82)
	#loc_40A = 33
	#loc_80A = 17
	print "loc_40A = %d, should be around 33"%loc_40A
	
	nf = len(ampList)
	print "nf = %d"%nf
	
	plt.cla()
	
	
	#plt.xlabel("Distance from origin in Fourier Space",fontsize=18)
	#plt.xlabel(r'Resolution ($1 /\AA$)',fontsize=18)
	plt.xlabel("Resolution",fontsize=18)
	plt.ylabel("Fourier Amplitude (Normalized)",fontsize=18)
	
	#colorList = ['g','r','g','c','grey','c','m','y']
	colorList = ['b','g','r','c','m','goldenrod','plum','lime','Violet','Hotpink','grey']
	#colorList = ['b','g','r','c','m','saddlebrown','goldenrod','plum','lime','Violet','Hotpink','grey']
	pp = []
	ll = []

	
	for i in range(nf):
		p, = plt.plot(xrange(len(ampList[i])),ampList[i],color=colorList[i],linewidth=1.0)
		peak40A_1st,peak40A_2nd,index40A_1st,index40A_2nd = find2peaks(loc_40A,ampList[i])
		peak80A_1st,peak80A_2nd,index80A_1st,index80A_2nd = find2peaks(loc_80A,ampList[i])
		# use plt.plot((x1, x2), (y1, y2), 'k-') to draw a line from the point (x1, y1) to the point (x2, y2) in color k
		#peak40A = (peak40A_1st+peak40A_2nd)/2
		#peak80A = (peak80A_1st+peak80A_2nd)/2
		peak40A = peak40A_1st
		peak80A = peak80A_1st
		if params['drawpeak']:
			#plt.plot((index40A_1st-5,index40A_1st+5),(peak40A,peak40A),linewidth=1.5,color=colorList[i])
			#plt.plot((index80A_1st-5,index80A_1st+5),(peak80A,peak80A),linewidth=1.5,color=colorList[i])
			order = 5
			v3 = max(ampList[i][index40A_1st*order-5:index40A_1st*order+5])
			p3 = ampList[i][index40A_1st*order-5:index40A_1st*order+5].index(v3)+index40A_1st*order-5
			plt.axvline(x=p3,color=colorList[i],linewidth=1.5)
			print "dimer spacing is %.2f"%(boxsize*apix*order/p3*2)
		pp.append(p)
		ll.append(params['fnames'].split(',')[i][:-4].replace("GTPgS",r'GTP$\gamma$S'))

		#if params['showres']:
			# FSC = 0.5
			#x_res1 = calcResNum(resList[i],fscList[i],0.5)
			#plt.axvline(x=x_res1,color=colorList[i],linestyle='dashed')
			#plt.figtext(0.15,0.8-0.1*i,"FSC = 0.5: %.2f A"%(1/x_res1),fontsize=14,color=colorList[i])
			# FSC = 0.143
			#x_res2 = calcResNum(resList[i],fscList[i],0.143)
			#plt.axvline(x=x_res2,color=colorList[i],linestyle='dashed')
			#plt.figtext(0.15,0.75-0.1*i,"FSC = 0.143: %.2f A"%(1/x_res2),fontsize=14,color=colorList[i])
			# Here 0.55 and 0.7 is position fraction
	
	# draw a horizontal line at y = 0.5
	#plt.axhline(y=0.5, color='k',linewidth=1.5, linestyle='dashed')
	#plt.axhline(y=0.143, color='k',linewidth=1.5, linestyle='dashed')
	plt.axhline(y=0,color='k',linewidth=1.5)

	
	# sets axes range
	res = params['res']
        rad = int(boxsize*apix/res)+2
	#rad = len(ampList[0])
	
	# set up the xticks
	#labels = [40, 40/2, 40/3, 40/4, 40/5, 40/6]
	labels = [40.0, 20.0, 13.3, 10.0, 8.0, 6.6, 5.7, 5.0, 4.4, 4.0]
	locs = [boxsize*apix/x for x in labels]
	plt.xticks(locs,labels,fontsize=16)
	plt.yticks(fontsize=16)

	x1,x2,y1,y2 = plt.axis()
        #plt.axis((0,rad,-1,16))
	plt.axis((0,rad,0,2.0))

	plt.legend(pp,ll,fontsize=18)
	#plt.title('Layer Line Intensities',fontsize=20)
	#if not os.path.exists('plot_LL.png'):
	#	plt.savefig("plot_LL.png")
	#else:
	#	plt.savefig("plot_LL_001.png")
	plt.savefig("plot_LL.png")
	#print "don't forget to rename plot.png"
	plt.show()

#===============
if  __name__ == "__main__":
	params = setupParserOptions()
	ampList = parseDataFiles(params)	
	plotData(ampList,params)

