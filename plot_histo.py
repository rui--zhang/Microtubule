#!/usr/bin/env python

import optparse
import itertools
import math
import os,sys
import glob
#import heapq
import numpy as np
import matplotlib.pyplot as plt
from numpy import arange

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <EB3clK_5.par>")
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="EB3clK_5.par")
	parser.add_option("-c",dest="column",type="int",metavar="INT", default=-1,
		help="choose which column to plot, 0-base")
        parser.add_option("-x",dest="xlabel",type="float",metavar="INT", default=14,
                help="plot a vertical line in x-axis")
	parser.add_option("--curve", action="store_true",dest="curve",default=False,
		help="plot as a curve, default plot as a histogram")
	parser.add_option("-t", dest="title", type="string", default='fig',
                help="give a title")

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

def mainloop(params):
	f1 = file(params['fpar'])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	column = params['column']
	data = [float(x.split()[column]) for x in l1]
	data = np.array(data)
	if params['curve']:
		n = len(data)
		xi = arange(0,n)-6
		plt.plot(xi, data, 'bo-')
		plt.xlabel("Protofilament Rotation",fontsize=18)
		plt.ylabel("Cross Correlation",fontsize=18)
	else:
		m = np.mean(data)
                std = np.std(data)
                plt.hist(data,50,range=[m-5*std,m+5*std],histtype='bar')
	xl = params['xlabel']
	plt.title('%s'%params['title'])
        plt.axvline(x=xl, color='k',linewidth=1.5, linestyle='dashed')
	plt.savefig("hist_%s.png"%params['fpar'][:-4])
	plt.show()

if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
