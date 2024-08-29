#!/usr/bin/env python
import optparse
import os,sys
#import numpy as np
import math
import itertools

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="13pf_1_r1.par")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT",default=1.32,
                help="pixel size in angstroms, default 1.32, no need to specify if fv9")
	parser.add_option("--fv",dest="fv",type="choice", metavar="['v8','v9']",
                choices=['v8','v9'],default='v9', help="input frealign format, default v9, output is always v9")
        parser.add_option("--greedy", action="store_true",dest="greedy",default=False,
                help="True if do frealign greedy")

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

def getMTlist(params):
	f1 = file(params["fpar"])
	l1 = f1.readlines()
	l2 = [x for x in l1 if x[0]!='C' and x!='\n']
	MTlist_dup = [abs(int(x.split()[7])) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	f1.close()
	return MTlist


# make phi follow a line
def unifyPhi(params,MT,philist,phi_MT):
	import copy
	MONOMER = 41.0
        PF = params['pf']
        RISE = MONOMER*3/PF
        TWIST = -360./PF
	if os.path.isfile('philist.txt'):
		#print "will use philist.txt\r",
		phi_MT,sign_MT = GetPhifromfile(MT,philist)
	else:
		#print "no philist.txt file specified, I will use median PHI"		
		#phi_MT = mymedian(philist)
		sign_MT = 1
	#print 'my phi is %.2f'%phi_MT
	n = len(philist)
	philist2 = copy.copy(philist)

	loc = philist.index(findclosest(phi_MT,philist))
	philist2[loc] = phi_MT
	#print "MT %d: loc = %.2f"%(MT,findclosest(phi_MT,philist))		

	for i in range(loc+1,n):
		dphi = philist2[i-1] - philist2[i]
		if abs(dphi) > 345.0:
			dphi = 0.0
		philist2[i] += TWIST*round(dphi/TWIST)

	for i in range(loc-1,-1,-1):
		dphi = philist2[i+1] - philist2[i]
		if abs(dphi) > 345.0:
			dphi = 0.0
		philist2[i] += TWIST*round(dphi/TWIST)

	return philist2,sign_MT

def sumAbs(mylist):
	mylist2 = [abs(i) for i in mylist]
	return sum(mylist2)
	
def shift40A_MT(shxlist,shylist,psi):
	MONOMER = 41.0
	shxlist_v1 = [i + MONOMER*math.cos(-psi/180.0*math.pi) for i in shxlist]
        shylist_v1 = [i + MONOMER*math.sin(-psi/180.0*math.pi) for i in shylist]
        shxlist_v2 = [i - MONOMER*math.cos(-psi/180.0*math.pi) for i in shxlist]
        shylist_v2 = [i - MONOMER*math.sin(-psi/180.0*math.pi) for i in shylist]

	if max(sumAbs(shxlist_v1),sumAbs(shylist_v1)) < max(sumAbs(shxlist_v2),sumAbs(shylist_v2)):
	#if max(abs(shx_v1),abs(shy_v1)) < max(abs(shx_v2),abs(shy_v2)):
                shxlist = shxlist_v1
                shylist = shylist_v1
        else:
                shxlist = shxlist_v2
                shylist = shylist_v2
	return shxlist,shylist

def getPhi_MT(l_MT):
	n = len(l_MT)
	
	psilist = [float(x.split()[1]) for x in l_MT]
	thetalist = [float(x.split()[2]) for x in l_MT]
	philist = [float(x.split()[3]) for x in l_MT]	
	shxlist = [float(x.split()[4]) for x in l_MT]
	shylist = [float(x.split()[5]) for x in l_MT]
	scorelist = [float(x.split()[14]) for x in l_MT]
	#MT = int(l_MT[0].split()[7])

        #start from the particle with the best score
        if params['fv'] == 'v8':
                loc = scorelist.index(min(scorelist))
        else:
                loc = scorelist.index(max(scorelist))

	phi_MT = philist[loc]
	psi_MT = psilist[loc]

	return phi_MT

def fv8tov9(l_MT,params):
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	l_MT_v9 = []
	apix = params['apix']
        for i in l_MT:
                t1 = i.split()
                count = int(t1[0])
                psi = float(t1[1])
                theta = float(t1[2])
                phi = float(t1[3])
                shx = float(t1[4])*apix
                shy = float(t1[5])*apix
		mag = float(t1[6])
                micro = int(t1[7])
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
                occ = 100.0
                logP = -22000.0
                sigma = 0.6
                score = 20.0
                change = 0.0
		l_MT_v9.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
	return l_MT_v9


def mainloop(params):
	f1 = file(params["fpar"])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	n = len(l1)
	MTlist = getMTlist(params)
	lastMT = int(l1[-1].split()[7])
	MTparlist = [[] for i in range(lastMT+1)]
	for i in l1:
		MT = int(i.split()[7])
		MTparlist[MT].append(i)

	if os.path.isfile('philist.txt'):
		print "will use philist.txt"
	else:
		print "no philist.txt file specified, I will use median PHI"
	
	fout2 = file("philist.txt","w")
	fout3 = file("MTlist.txt","w")
	MTnew = 0
	for MT in MTlist:
		print "working on MT %d\t\r"%MT,
		l_MT = MTparlist[MT]
		if params['greedy']:
			fout = file("MT_%d.par"%MT,"w")
			fout.writelines(l_MT)
			fout.close()
		if params['fv'] == 'v8':
			l_MT = fv8tov9(l_MT,params)
		NN = len(l_MT)
		#if NN > 300:
		#	fout.writelines(unifyAll(params,l_MT[:NN/2],MT))
		#	fout.writelines(unifyAll(params,l_MT[NN/2:],MT))
		#else:
		#	fout.writelines(unifyAll(params,l_MT,MT))
		phi_MT = getPhi_MT(l_MT)
		fout2.write("%d\t%.2f\n"%(MT,phi_MT))
		fout3.write("%d\n"%(MT))
			
	f1.close()
	fout2.close()
	fout3.close()

if __name__ == "__main__":
	params = setupParserOptions()
	os.system("rm -f philist.txt")
	mainloop(params)
