#!/usr/bin/env python
import optparse
import os,sys
import numpy as np
import math
import itertools
from random import randint
import time
import os.path

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="13pf_1_r1.par")
	parser.add_option("--pf", dest="pf", type="int", metavar="INT",
                help="13 or 14")
	parser.add_option("--iterr", dest="iterr", type="int", metavar="INT",
                help="current iteration")
	#parser.add_option("--uniPhi",dest="uniPhi",action="store_true",default=True,
        #        help="if True, unify Phi; default=False")


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

def closer2origin(shx,shy,psi):
	DIMER = 82.0
	shx_new1 = shx + DIMER*math.cos(-psi/180.0*math.pi)
	shy_new1 = shy + DIMER*math.sin(-psi/180.0*math.pi)		
	shx_new2 = shx - DIMER*math.cos(-psi/180.0*math.pi)
	shy_new2 = shy - DIMER*math.sin(-psi/180.0*math.pi)
	if max(abs(shx_new1),abs(shy_new1)) < max(abs(shx),abs(shy)):
		shx = shx_new1
		shy = shy_new1
	if max(abs(shx_new2),abs(shy_new2)) < max(abs(shx),abs(shy)):
		shx = shx_new2
		shy = shy_new2
	return shx,shy


def mymedian(mylist):
	import copy
	n = len(mylist)
	if n%2 == 1:
		return np.median(mylist)
	else:
		mylist2 = copy.copy(mylist)
		mylist2.sort()
		return np.median(mylist2[1:])

def getMTlist(params):
	f1 = file(params["fpar"])
	l1 = f1.readlines()
	l2 = [x for x in l1 if x[0]!='C' and x!='\n']
	MTlist_dup = [abs(int(x.split()[7])) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	f1.close()
	return MTlist

def closer2origin(shx,shy,psi):
        DIMER = 82.0
        shx_new1 = shx + DIMER*math.cos(-psi/180.0*math.pi)
        shy_new1 = shy + DIMER*math.sin(-psi/180.0*math.pi)
        shx_new2 = shx - DIMER*math.cos(-psi/180.0*math.pi)
        shy_new2 = shy - DIMER*math.sin(-psi/180.0*math.pi)
        if max(abs(shx_new1),abs(shy_new1)) < max(abs(shx),abs(shy)):
                shx = shx_new1
                shy = shy_new1
        if max(abs(shx_new2),abs(shy_new2)) < max(abs(shx),abs(shy)):
                shx = shx_new2
                shy = shy_new2
        return shx,shy


def freali_exe(path,line,pf,iterr):
	os.chdir(path)
	os.system("cp ../%spf_%d_r1.par %spf_%d_r1.par_orig"%(pf,iterr-1,pf,iterr-1))
	f2 = file('%spf_%d_r1.par_orig'%(pf,iterr-1))
	l2 = f2.readlines()
	#l2 = [x for x in ll2 if x[0]!='C' and x!='\n']
	f2.close()

	start = int(line.split()[0])
	l2_new = l2[:start-1]+[line]+l2[start:]
	#print start
	#print l2_new[start-2:start+3]

	fout = file('%spf_%d_r1.par'%(pf,iterr-1),"w")
	fout.writelines(l2_new)

        fout.close()

	os.system('cp ../mparameters .')
	os.system('cp ../mult_hrefine_n.com .')
	os.system('ln -sf ../%spf_%d_r1.mrc .'%(pf,iterr-1))

	start = int(line.split()[0])
	print "./mult_hrefine_n.com %d %d %d 1"%(start,start,iterr)
	os.system("./mult_hrefine_n.com %d %d %d 1"%(start,start,iterr))
	f1 = file('%spf_%d_r1.par_%d_%d'%(pf,iterr,start,start))
	l1 = f1.readlines()
	f1.close()
	print l1[-3]
	os.chdir("..")
	return l1[-3]

def copy_params(line_new,line_old):
	#print 'new line is %s'%line_new
	#print 'old line is %s'%line_old
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	t1 = line_new.split()
	t2 = line_old.split()
	count = int(t2[0])
        psi = float(t1[1])
        theta = float(t1[2])
        phi = float(t1[3])
        shx = float(t1[4])
        shy = float(t1[5])
	#shx,shy = closer2origin(shx,shy,psi)
	#shx,shy = closer2origin(shx,shy,psi)
        mag = float(t2[6])
        micro = int(t2[7])
        df1 = float(t2[8])
        df2 = float(t2[9])
        angast = float(t2[10])
	occ = float(t2[11])
       	logP = float(t2[12])
        sigma = float(t2[13])
        score = float(t2[14])
        change = float(t2[15])
	if micro > 99999:
		FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d %6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	return FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change)

def unifyAll(path,params,l_MT,MT):
	thresh = 1.5
	pf = params['pf']
	iterr = params['iterr']
	n = len(l_MT)
	scorelist = [float(x.split()[14]) for x in l_MT]
	philist = [float(x.split()[3]) for x in l_MT]
	psilist = [float(x.split()[1]) for x in l_MT]
	#MT = int(l_MT[0].split()[7])

	l_MT_new = [[] for x in l_MT]
	l_MT_new_tmp = [[] for x in l_MT]
	l_MT_new_tmp2 = [[] for x in l_MT]
		
	#start from the particle with the best score
	loc = scorelist.index(max(scorelist))
	
	bestscore = scorelist[loc]
	bestphi = philist[loc]
	bestpsi = psilist[loc]
	#print "bestscore is %.2f"%bestscore
	
	l_MT_new[loc] = l_MT[loc]
	
	# process from loc to the end
	for i in range(loc+1,n):
		#if abs(scorelist[i]-bestscore) < 1.5 or (scorelist[i]>bestscore):
		if scorelist[i] > bestscore-thresh and abs(bestphi-philist[i]) < 5.0:
			l_MT_new[i] = l_MT[i]
		else:
			l_MT_new_tmp[i] = copy_params(l_MT_new[i-1],l_MT[i])
			l_MT_new_tmp2[i] = freali_exe(path,l_MT_new_tmp[i],pf,iterr)
			bestscore_tmp = float(l_MT_new_tmp2[i].strip().split()[-2])
			# only update if the score is better
			if bestscore_tmp > scorelist[i]:
				l_MT_new[i] = l_MT_new_tmp2[i]
                        elif abs(bestpsi-psilist[i]) > 5.0:
                                l_MT_new[i] = l_MT_new_tmp2[i]
			else:
				l_MT_new[i] = l_MT[i]
		bestscore = float(l_MT_new[i].strip().split()[-2])
		bestphi = float(l_MT_new[i].strip().split()[3])

	# reset bestscore and bestphi
	bestscore = scorelist[loc]
	bestphi = philist[loc]
	bestpsi = psilist[loc]

	for i in range(loc-1,-1,-1):
		#if abs(scorelist[i]-bestscore) < 1.5 or (scorelist[i]>bestscore):
		if scorelist[i] > bestscore-thresh and abs(bestphi-philist[i]) < 5.0:
			l_MT_new[i] = l_MT[i]
		else:
			l_MT_new_tmp[i] = copy_params(l_MT_new[i+1],l_MT[i])
			l_MT_new_tmp2[i] = freali_exe(path,l_MT_new_tmp[i],pf,iterr)
			bestscore_tmp = float(l_MT_new_tmp2[i].strip().split()[-2])
			# only update if the score is better
			if bestscore_tmp > scorelist[i]:
				l_MT_new[i] = l_MT_new_tmp2[i]
			elif abs(bestpsi-psilist[i]) > 5.0:
				l_MT_new[i] = l_MT_new_tmp2[i]
			else:
				l_MT_new[i] = l_MT[i]
		bestscore = float(l_MT_new[i].strip().split()[-2])
		bestphi = float(l_MT_new[i].strip().split()[3])

	return l_MT_new


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

	fout = file("%s-greedy"%params["fpar"],"w")
	fout2 = file('finished.txt',"a")
	
	for MT in MTlist:
		print "working on MT %d\t\r"%MT,
		path = "zdir_MT_%d"%MT
		if not os.path.exists(path):
    			os.makedirs(path)
		l_MT = MTparlist[MT]
		l_MT_new2 = unifyAll(path,params,l_MT,MT)
		fout.writelines(l_MT_new2)
		os.system("mv MT_%d.par* finished/"%(MT))
		os.system("rm -rf zdir_MT_%d"%MT)
		if len(l_MT) == len(l_MT_new2): 
			fout2.write("%s\n"%MT)
			
	f1.close()
	fout.close()
	fout2.close()

if __name__ == "__main__":
	params = setupParserOptions()
	if not params['fpar']:
                sys.exit()
	mainloop(params)
