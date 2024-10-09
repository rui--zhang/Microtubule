#!/usr/bin/env python
import optparse
import os,sys
import numpy as np
import math
import itertools


def setupParserOptions():
	parser = optparse.OptionParser()
	#parser.set_usage("%prog -f <4TuIntra_pf.pdb>")
	parser.add_option("-f",dest="fpdb",type="string",metavar="FILE", default='4TuIntra_pf.pdb',
		help="4TuIntra_pf.pdb")
	
	
	options,args = parser.parse_args()

	if len(args) > 1:
		parser.error("Unknown commandline options: " +str(args))

	#if len(sys.argv) < 2:
	#	parser.print_help()
	#	sys.exit()

	params={}

	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params


def dist(A,B):
	return math.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2+(A[2]-B[2])**2)
	

def diff(m,n):
	zm = m[2]
	zn = n[2]
	return int(zm-zn)
	

def dimer(params):
	fpdb = params['fpdb']
	fout = file('dist.txt',"w")
	f1 = file(fpdb)
        ll1 = f1.readlines()
        l1 = [i for i in ll1 if (i[:4]=='ATOM')]

	lB = [x for x in l1 if x.split()[4]=='B' and x.split()[2]=='CA']
	lD = [x for x in l1 if x.split()[4]=='D' and x.split()[2]=='CA']
	cordlistB = []
	cordlistD = []
	distlist = []

	for i in lB:
		tmp = i.split()
		resi = int(tmp[5])
		cordlistB.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])
	
	for i in lD:
		tmp = i.split()
		resi = int(tmp[5])
		cordlistD.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])

	#print len(cordlistA)
	#print len(cordlistB)
	if len(cordlistB) != len(cordlistD):
		print "len(cordlistB) != len(cordlistD)"
		sys.exit()
	
	for i in range(len(cordlistB)):
		dd = dist(cordlistB[i],cordlistD[i])
		distlist.append(dd)
		fout.write("%d\t%.2f\n"%(i,dd))
	
	distlist.sort()

	print "dimer"
	#print "median"
	print np.median(distlist)
	#print "average"
	#print np.average(distlist)
	#print distlist

	f1.close()
	fout.close()


def intra(params):
	fpdb = params['fpdb']
	fout = file('dist.txt',"w")
	f1 = file(fpdb)
        ll1 = f1.readlines()
        l1 = [i for i in ll1 if (i[:4]=='ATOM')]

	lA = [x for x in l1 if x.split()[4]=='A' and x.split()[2]=='CA']
	lB = [x for x in l1 if x.split()[4]=='B' and x.split()[2]=='CA']
	cordlistA = []
	cordlistB = []
	distlist = []

	for i in lA:
		tmp = i.split()
		resi = int(tmp[5])
		if resi < 38 or (resi > 47 and resi < 171):
		#if resi < 38 or (resi > 47 and resi < 361) or (resi > 368 and resi < 437):
			cordlistA.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])
	
	for i in lB:
		tmp = i.split()
		resi = int(tmp[5])
		if resi < 38 or (resi > 45 and resi < 169):
		#if resi < 38 or (resi > 45 and resi < 359) or (resi > 358 and resi < 427):
			cordlistB.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])

	#print len(cordlistA)
	#print len(cordlistB)
	if len(cordlistA) != len(cordlistB):
		print "len(cordlistA) != len(cordlistB)"
		sys.exit()
	
	for i in range(len(cordlistA)):
		dd = dist(cordlistA[i],cordlistB[i])
		distlist.append(dd)
		fout.write("%d\t%.2f\n"%(i,dd))
	
	distlist.sort()

	print "intra"
	#print "median"
	print np.median(distlist)
	#print "average"
	#print np.average(distlist)
	#print distlist

	f1.close()
	fout.close()



def inter(params):
	fpdb = params['fpdb']
	fout = file('dist.txt',"w")
	f1 = file(fpdb)
        ll1 = f1.readlines()
        l1 = [i for i in ll1 if (i[:4]=='ATOM')]

	lC = [x for x in l1 if x.split()[4]=='C' and x.split()[2]=='CA']
	lB = [x for x in l1 if x.split()[4]=='B' and x.split()[2]=='CA']
	cordlistC = []
	cordlistB = []
	distlist = []


	for i in lC:
		tmp = i.split()
		resi = int(tmp[5])
		if resi < 38 or (resi > 47 and resi < 171):
		#if resi < 38 or (resi > 47 and resi < 241):
		#if resi < 38 or (resi > 47 and resi < 361) or (resi > 368 and resi < 437):
			cordlistC.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])
	
	for i in lB:
		tmp = i.split()
		resi = int(tmp[5])
		if resi < 38 or (resi > 45 and resi < 169):
		#if resi < 38 or (resi > 45 and resi < 239):
		#if resi < 38 or (resi > 45 and resi < 359) or (resi > 358 and resi < 427):
			cordlistB.append([float(tmp[-6]),float(tmp[-5]),float(tmp[-4])])

	if len(cordlistC) != len(cordlistB):
		print "len(cordlistC) != len(cordlistB)"
		sys.exit()
	
	for i in range(len(cordlistC)):
		dd = dist(cordlistC[i],cordlistB[i])
		distlist.append(dd)
		fout.write("%d\t%.2f\n"%(i,dd))
	
	distlist.sort()

	print "inter"	
	print np.median(distlist)
	#print distlist

	f1.close()
	fout.close()



if __name__ == "__main__":
	params = setupParserOptions()
	intra(params)
	inter(params)
	dimer(params)
