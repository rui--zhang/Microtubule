#!/usr/bin/env python

import os,sys
import optparse
from EMAN2 import *
import math
#from scipy import stats
import numpy as np
from numpy import arange
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import matplotlib.pyplot as plt


#===============
def setupParserOptions():
	parser = optparse.OptionParser()
	
	#parser.set_usage("%prog -f fscData1.txt,fscData2.txt -c 3")
	parser.add_option("-f", dest="fname", type="string", default='perptcl_stats.txt',
		help="default perptcl_stats.txt")
	parser.add_option("--noplotfig", action="store_true",dest="noplotfig",default=False,
                help="plot figure, for debugging")
        parser.add_option("--label",dest="label",type="string",
                help="label")
	
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

def getMTlist(l1):
	MTlist_dup = [int(x.split()[2]) for x in l1 if x!='\n']
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	return MTlist

def parse1MT(lines_MT):
	n = len(lines_MT)/13
	pres1list = np.array([0.0]*13)
	pres2list = np.array([0.0]*13)
	dpreslist = np.array([0.0]*13)
	for i in range(n):
		pres1 = np.array([float(x.split()[-2]) for x in lines_MT[i*13:(i+1)*13]])
		pres2 = np.array([float(x.split()[-1]) for x in lines_MT[i*13:(i+1)*13]])
		dpres = np.array([float(x.split()[-3]) for x in lines_MT[i*13:(i+1)*13]])
		pres1list += pres1
		pres2list += pres2
		dpreslist += dpres
	return pres1list/n,pres2list/n,dpreslist/n

def CutbyPercentage(mylist,perc):
        import copy
        n = len(mylist)
        mylist2 = copy.copy(mylist)
        mylist2.sort()
        cutoff = mylist2[int(perc*n)-1]
        return cutoff

def plot1ptcl(pres1list,pres2list,dpreslist,MT,params):
	plt.ioff()
	plt.close('all')
	#pres1 = [float(x.split()[-2]) for x in l1]
	#pres2 = [float(x.split()[-1]) for x in l1]
	#dpres = [float(x.split()[-3]) for x in l1]

	xi = arange(-6,7)
	pres2list_fit = np.polyfit(xi,pres2list,2)
	pres1list_fit = np.polyfit(xi,pres1list,2)
	dpreslist_fit = np.polyfit(xi,dpreslist,2)
	
	peak = -dpreslist_fit[1]/2/dpreslist_fit[0]
	try:
		val = int(peak)
	except:
		peak = -999
	## take the max(abs)
	if peak < -6 or peak >6:
		if abs(np.polyval(dpreslist_fit,-6)) >= abs(np.polyval(dpreslist_fit,6)):
			MTpeak = -6
		else:
			MTpeak = 6
	else:
		MTpeak = int(round(peak))
	peakdpres = np.polyval(dpreslist_fit,MTpeak)
	
	if not params['noplotfig']:
		#MT = int(l1[0].split()[2])
		plt.subplot(221)
		plt.plot(xi,pres1list,color='b',linewidth=2)
		plt.plot(xi,pres2list,color='r',linewidth=2)
		plt.plot(xi,np.polyval(pres1list_fit,xi),'b-')
		plt.plot(xi,np.polyval(pres2list_fit,xi),'r-')
		plt.xlim([-6,6])
		plt.ylabel("Similarity Score",fontsize=18)
		plt.xlabel('Candidate Seam Location',fontsize=18)
		plt.tick_params(axis='both', which='major', labelsize=18)
		plt.subplot(223)
		plt.plot(xi,dpreslist,color='g',linewidth=2)
		# This is not needed, dpres_fit2 should exactly match dpres_fit
		#dpres_fit2 = np.polyval(pres1_fit,xi)-np.polyval(pres2_fit,xi)
		#plt.plot(xi,dpres_fit2,'m-')
		plt.plot(xi,np.polyval(dpreslist_fit,xi),'g-')
		# for y = ax2+bx+c, peak = -b/(2a)
		peak = -dpreslist_fit[1]/2/dpreslist_fit[0]
		plt.axvline(x=peak,color='g',linestyle='dashed')
		plt.ylabel("Similarity Score",fontsize=18)
		plt.xlabel('Candidate Seam Location',fontsize=18)
		#plt.axvline(x=-4,color='k',linestyle='dashed')
		plt.axhline(y=0,color='k',linestyle='dashed')
		plt.figtext(0.12,0.9,'Score1',fontsize=16,color='b')
		plt.figtext(0.42,0.9,'Score2',fontsize=16,color='r')
		plt.figtext(0.2,0.2,'Score1-Score2',fontsize=16,color='g')
		plt.xlim([-6,6])
		#plt.ylim([0,10])
		plt.tick_params(axis='both', which='major', labelsize=18)

		plt.tight_layout()
		#plt.legend(ll)
		plt.savefig("plot_seam/MT%d_%s.png"%(MT,params['label']),dpi=50)
		#plt.show()

	return MTpeak,peakdpres

def mainloop(params):
	if not os.path.exists('plot_seam'):
		os.makedirs('plot_seam')
	f1 = file(params['fname'])
	l1 = f1.readlines()
	MTlist = getMTlist(l1)
	lastMT = MTlist[-1]
	MTparlist = [[] for i in range(lastMT+1)]

	for i in l1:
		if i == '\n':
			continue
		MT = int(i.split()[2])
		MTparlist[MT].append(i)

	fout = file('picklist3.txt',"w")

	peakdpresList = []

        from progressbar import ProgressBar
        pbar = ProgressBar()

	for MT in pbar(MTlist):
		pres1list,pres2list,dpreslist = parse1MT(MTparlist[MT])
		dpreslist_abs = abs(dpreslist)
		MTpeak, peakdpres = plot1ptcl(pres1list,pres2list,dpreslist,MT,params)
		MTpeak = list(dpreslist_abs).index(max(dpreslist_abs))
		peakdpres = dpreslist[MTpeak]
		
		fout.write("%d\t%d\t%.2f\n"%(MT,MTpeak-6,peakdpres))
		peakdpresList.append(abs(peakdpres))

	cutoff = CutbyPercentage(peakdpresList,0.2)
	print "discard 20 percent of MTs, dpres cutoff = %.2f"%cutoff
	fout2 = file('goodMTlist_cut%.1f.txt'%cutoff,"w")
	n = len(MTlist)
	for i in range(n):
		MT = MTlist[i]
		peakdpres = peakdpresList[i]
		if abs(peakdpres) > cutoff:	
			fout2.write("%d\t%.2f\n"%(MT,peakdpres))	

	f1.close()
	fout.close()
	fout2.close()

#===============
if  __name__ == "__main__":
	params = setupParserOptions()
	os.system("rm -rf plot_seam")
	mainloop(params)

