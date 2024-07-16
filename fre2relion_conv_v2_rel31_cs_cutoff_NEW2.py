#!/usr/bin/env python

import optparse
import os,sys
#from EMAN2 import *
import math


def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <EB3clK_5.par>")
	parser.add_option("--f1",dest="fpar1",type="string",metavar="FILE",
		help="frealign 14pf_2_r1.par")
	parser.add_option("--f2",dest="fpar2",type="string",metavar="FILE",
		help="relion particles_reorder.star, from Particle Extraction, and after re-order")
	parser.add_option("--cut", dest="cut", type="float", metavar="FLOAT",
                help="cutoff value")

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

def getColNum(lines,word):
	for i in lines:
		if word in i:
			return int(i.split()[-1][1:])

def fre2relion(params):
	f1 = file(params['fpar1'])
	ll1 = f1.readlines()
	l1 = [x for x in ll1 if x[0]!='C' and x!='\n']

	f2 = file(params['fpar2'])
	ll2 = f2.readlines()
	l2 = [x for x in ll2 if len(x.split())>11]

	cutoff = params['cut']

	fout = file('%s_fre2relion_above%d_NEW2.star'%((params['fpar2'][:-5]),cutoff),"w")
	fout2 = file('%s_fre2relion_below%d_NEW2.star'%((params['fpar2'][:-5]),cutoff),"w")
	

	n1 = len(l1)
	n2 = len(l2)

	print n1,n2
	if n1 != n2:
		print "Error, n1 != n2"

	# first write out the header
	for i in range(100):
		t2 = ll2[i].split()
		if len(t2) < 15:
			fout.write("%s"%ll2[i])
			fout2.write("%s"%ll2[i])
	
	#if (num_df1 != 10) or (num_df2 != 11) or (num_angast != 12) or (num_amp != 17):
	#	print "Error, columns don't match!"

#_rlnImageName #1
#_rlnMicrographName #2
#_rlnCoordinateX #3
#_rlnCoordinateY #4
#_rlnAngleRot #5
#_rlnAngleTilt #6
#_rlnAnglePsi #7
#_rlnOriginXAngst #8
#_rlnOriginYAngst #9
#_rlnDefocusU #10
#_rlnDefocusV #11
#_rlnDefocusAngle #12
#_rlnPhaseShift #13
#_rlnCtfBfactor #14
#_rlnOpticsGroup #15
#_rlnRandomSubset #16
#_rlnClassNumber #17
#_rlnHelicalTubeID #18

	for i in range(n2):
		t1 = l1[i].split()
		psi = float(t1[1])
		theta = float(t1[2])
		phi = float(t1[3])
		shx = float(t1[4])
		shy = float(t1[5])
		df1 = float(t1[8])
		df2 = float(t1[9])
		angast = float(t1[10])
		pres = float(t1[-2])
		t2 = l2[i].split()
#		fout.write('%s\t%s\t%s\t%s\t%s\t %.6f\t%.6f\t%.6f\t %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n'%(t2[0],t2[1],t2[2],t2[3],t2[4],df1,df2,angast,t2[8],t2[9],t2[10],t2[11],t2[12],t2[13],t2[14],t2[15],t2[16],1,phi,theta,psi,-1*shx/apix,-1*shy/apix))
		#fout.write('%s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\t%d\t%.6f\t %.6f\t%.6f\t%.6f\t%.6f\n'%(t2[0],t2[1],t2[2],t2[3],t2[4],t2[5],t2[6],t2[7],t2[8],t2[9],t2[10],t2[11],t2[12],1,phi,theta,psi,-1*shx/apix,-1*shy/apix))
		if pres >= cutoff:
			fout.write('%s\t%s\t%s\t%s\t%.6f\t %.6f\t%.6f\t%.6f\t%.6f\t%s\t %s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\n'%(t2[0],t2[1],t2[2],t2[3],phi,theta,psi,-1*shx,-1*shy,t2[9],t2[10],t2[11],t2[12],t2[13],t2[14],t2[15],t2[16],t2[17]))
		else:
			fout2.write('%s\t%s\t%s\t%s\t%.6f\t %.6f\t%.6f\t%.6f\t%.6f\t%s\t %s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\n'%(t2[0],t2[1],t2[2],t2[3],phi,theta,psi,-1*shx,-1*shy,t2[9],t2[10],t2[11],t2[12],t2[13],t2[14],t2[15],t2[16],t2[17]))


	f1.close()
	f2.close()
	fout.close()
	fout2.close()

	os.system("wc %s_fre2relion_above%d_NEW2.star"%((params['fpar2'][:-5]),cutoff))
	os.system("wc %s_fre2relion_below%d_NEW2.star"%((params['fpar2'][:-5]),cutoff))

if __name__ == "__main__":
	params = setupParserOptions()
	fre2relion(params)
