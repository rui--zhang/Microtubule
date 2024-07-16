#!/usr/bin/env python
from EMAN2 import *
import math
import optparse
import os,sys
import glob
import numpy as np

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -f <EB3clK_5.par>")
        parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
                help="EB3clK_5.par")
	parser.add_option("-s",dest="stack",type="string",metavar="FILE",
		help="start.hed")
	#parser.add_option("--pf", dest="pf", type="int", metavar="int",
        #        help="PF")
	parser.add_option("--outmrc", action="store_true",dest="outmrc",default=False,
		help="output to mrc format")
	parser.add_option("--offsetshxy", action="store_true",dest="offsetshxy",default=False,
		help="offset shxy")
	parser.add_option("--apix", dest="apix", type="float", metavar="FLOAT",
                help="pixel size in angstroms")

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

def getMTlist(l1):
	l2 = [x for x in l1 if x[0]!='C' and x!='\n']
	MTlist_dup = [int(x.split()[7]) for x in l2]
	MTlist = list(set(MTlist_dup))
	MTlist.sort()
	return MTlist

def grepMT(l1,MT):
	#f1 = file(fname)
	#l1 = f1.readlines()
        l2 = [x for x in l1 if x[0]!='C' and x!='\n']
        l3 = [x for x in l2 if int(x.split()[7])==MT]
	#f1.close()
        return l3

def mymedian(mylist):
	import copy
	n = len(mylist)
	if n%2 == 1:
		return np.median(mylist)
	else:
		mylist2 = copy.copy(mylist)
		mylist2.sort()
		return np.median(mylist2[1:])

def convertSameLine(l_MT,params):
	apix = params['apix']
	if params['pf'] == 12:
		RISE = 10.44/apix
		TWIST = -29.86
	elif params['pf'] == 13:
		RISE = 9.5/apix
		TWIST = -27.71
	elif params['pf'] == 14:
		RISE = 8.85/apix
		TWIST = -25.76
	else:
		print "please specify PF"
		sys.exit()
	MONOMER = 41.0/apix
	
	n = len(l_MT)
	
	psilist = [float(x.split()[1]) for x in l_MT]
	thetalist = [float(x.split()[2]) for x in l_MT]
	philist = [float(x.split()[3]) for x in l_MT]	
	shxlist = [float(x.split()[4]) for x in l_MT]
	shylist = [float(x.split()[5]) for x in l_MT]
	MT = int(l_MT[0].split()[7])

	#if os.path.isfile('picklist.txt'):
	#	print "will use picklist.txt\r",
	#	phi_MT,sign_MT = GetPhi(MT,philist)
	#else:
	#	phi_MT = np.median(philist)
	#	sign_MT = 1
	
	for i in range(n):
		dphi = phi_MT - philist[i]
		phi_new = philist[i] + TWIST*round(dphi/TWIST)
		if abs(phi_new-philist[i]) > 15.0:
			#scale = (phi_new-philist[i])/TWIST
			scale = round(dphi/TWIST)
			shxlist[i] += scale*RISE*math.cos(-psilist[i]/180.0*math.pi)
			shylist[i] += scale*RISE*math.sin(-psilist[i]/180.0*math.pi)
			philist[i] = phi_new
		# shift 1 monomer
		if sign_MT < 0:
			shxlist[i] += MONOMER*math.cos(-psilist[i]/180.0*math.pi)
			shylist[i] += MONOMER*math.sin(-psilist[i]/180.0*math.pi)

	loc = findMinShxy(shxlist,shylist)
	#print "MT = %d, loc = %d"%(MT,loc)
	for i in range(loc+1,n):
		shxlist[i],shylist[i] = findClosestShxy(shxlist[i],shylist[i],shxlist[i-1],shylist[i-1],psilist[i],apix)
	for i in range(loc-1,-1,-1):
		shxlist[i],shylist[i] = findClosestShxy(shxlist[i],shylist[i],shxlist[i+1],shylist[i+1],psilist[i],apix)
	
	l_MT_new = []
	#FORMAT = "%7d%8.2f%8.2f%8.2f%8.2f%8.2f%7.0f.%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n"
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	for i in range(n):
		t1 = l_MT[i].split()
		count = int(t1[0])
                #psi = float(t1[1])
                #theta = float(t1[2])
                #phi = float(t1[3])
                #shx = float(t1[4])
                #shy = float(t1[5])
		psi = psilist[i]
		theta = thetalist[i]
		phi = philist[i]
		shx = shxlist[i]
		shy = shylist[i]
                mag = float(t1[6])
                micro = int(t1[7])
                df1 = float(t1[8])
                df2 = float(t1[9])
                angast = float(t1[10])
                occ = float(t1[11])
                logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
                l_MT_new.append(FORMAT%(count,psi,theta,phi,shx,shy,mag,micro,df1,df2,angast,occ,logP,sigma,score,change))
	return l_MT_new

def getoffset(mylist):
	aboffset = 999
	offset = 999
	for i in mylist:
		if abs(i) < aboffset:
			aboffset = abs(i)
			offset = i
	return offset

def MakeMTavg(l_MT,stack,nx,offsetshxy,apix):
	start = int(l_MT[0].split()[0])
	end = int(l_MT[-1].split()[0])
	MTstack = EMData.read_images(stack,range(start-1,end))
	
	MTavg = EMData(nx,nx)
	MTavg.to_zero()
	n = len(l_MT)
	
	philist = []
	thetalist = []
	logPlist = []
	sigmalist = []
	scorelist = []
	changelist = []
	
	
	if offsetshxy:
		shxlist = [float(x.split()[4]) for x in l_MT]
		shylist = [float(x.split()[5]) for x in l_MT]
		shx_offset = getoffset(shxlist)
		shy_offset = getoffset(shylist)
	else:
		shx_offset = 0
		shy_offset = 0
	
	for i in range(n):
		img = MTstack[i]
		t1 = l_MT[i].split()
		psi = float(t1[1])
		theta = float(t1[2])
		phi = float(t1[3])
		shx = float(t1[4])-shx_offset
		shy = float(t1[5])-shy_offset
		mag = float(t1[6])
		micro = int(t1[7])
		df1 = float(t1[8])
		df2 = float(t1[9])
		angast = float(t1[10])
                #occ = float(t1[11])
                logP = float(t1[12])
                sigma = float(t1[13])
                score = float(t1[14])
                change = float(t1[15])
		#
		philist.append(phi)
		thetalist.append(theta)
		logPlist.append(logP)
		sigmalist.append(sigma)
		scorelist.append(score)
		changelist.append(change)
		# Frealign applies the shifts first and then the rotations.
		t1 = Transform()
		t1.set_trans(-shx/apix,-shy/apix)
		img2 = img.process("xform",{"transform":t1})
		t2 = Transform()
		t2.set_rotation({"type":"2d","alpha":-psi})
		img3 = img2.process("xform",{"transform":t2})
		MTavg.add(img3)
	
	phi_MT = mymedian(philist)
	theta_MT = mymedian(thetalist)
	logP_MT = np.median(logPlist)
	sigma_MT = np.median(sigmalist)
	score_MT = np.median(scorelist)
	change_MT = np.median(changelist)
	#FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	FORMAT2 = "%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	l_MTavg = FORMAT2%(0,theta_MT,phi_MT,0,0,mag,micro,df1,df2,angast,100,logP_MT,sigma_MT,score_MT,change_MT)
	return MTavg,l_MTavg

def mainloop(params):
	apix = params['apix']
	f1 = file(params["fpar"])
        ll1 = f1.readlines()
        l1 = [x for x in ll1 if x[0]!='C' and x!='\n']
	
	MTlist = getMTlist(l1)
	#print MTlist
	#nMT = len(MTlist)

	lastMT = int(l1[-1].split()[7])
	#print lastMT
	MTparlist = [[] for i in range(lastMT+1)]
	for i in l1:
		MT = int(i.split()[7])
		MTparlist[MT].append(i)
	
	# first calculate the num of particles
	count = 0
	SEG = 7
	for MT in MTlist:
		#print "pre-working on MT %d\t\r"%MT,
		#l_MT = grepMT(l1,MT)
		l_MT = MTparlist[MT]
		nptcl = len(l_MT)
		if nptcl < 3:
			continue
		
		nSEG = nptcl/SEG
		if nptcl%SEG > SEG-2 or nSEG == 0:
			nSEG += 1
		count += nSEG
	
	stack = params["stack"]
	im = EMData(stack,0)
	nx = im.get_xsize()
	del im
	print "# of MTs: %d"%count
	if params['outmrc']:
		print "Allocating space for MTavgstack...\n"
		MTavgstack = EMData(nx,nx,count)
		MTavgstack.write_image("MTavgstack_break.mrc")
		print "Done allocating"

	f2 = file("%s_MTavg_break"%params["fpar"],"w")
	#FORMAT = "%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.0f%6d%9.1f%9.1f%8.2f%7.2f%8.2f\n"
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	
	if params['offsetshxy']:
		offsetshxy = 1
	else:
		offsetshxy = 0

	# reset count = 0
	count = 0

        from progressbar import ProgressBar
        pbar = ProgressBar()

	for MT in pbar(MTlist):
		#print "working on MT %d\t\r"%MT,
		#l_MT = grepMT(l1,MT)
		l_MT = MTparlist[MT]
		nptcl = len(l_MT)
		if nptcl < 3:
			continue
		tmp0 = l_MT[0].split()
		mag = float(tmp0[6])
		#df1 = float(tmp0[8])
		#df2 = float(tmp0[9])
		angast = float(tmp0[10])
		pres = float(tmp0[11])
		
		nSEG = nptcl/SEG
		ttt2 = nptcl%SEG
		if nptcl%SEG > SEG-2 or nSEG == 0:
			nSEG += 1
		for i in range(nSEG):
			start = i*SEG
			#end = (i+1)*SEG
			if i == nSEG-1:
				end = nptcl
			else:
				end = (i+1)*SEG
			l_MT_SEG = l_MT[start:end]
			#l_MT_SEG_new = convertSameLine(l_MT_SEG,params)
			MTavg,l_MTavg = MakeMTavg(l_MT_SEG,stack,nx,offsetshxy,apix)
			if params['outmrc']:		    
				region = Region(0,0,count,nx,nx,1)
				MTavg.write_image("MTavgstack_break.mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
			else:
				MTavg.write_image("MTavgstack_break.hed",-1)
			f2.write("%7d%s"%(count+1,l_MTavg))
			count += 1

	f1.close()
	f2.close()


if __name__ == "__main__":
        params = setupParserOptions()
	if params['outmrc']:
		try:
			os.remove("MTavgstack_break.mrc")
		except:
			pass
	else:
		try:
                        os.remove("MTavgstack_break.hed")
			os.remove("MTavgstack_break.img")
                except:
                        pass
	#for i in glob.glob("MTavgstack_break.*"):
        #        os.remove(i)
	mainloop(params)
