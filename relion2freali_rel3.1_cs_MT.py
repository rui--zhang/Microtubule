#!/usr/bin/env python

import optparse
import os,sys
#from EMAN2 import *
import math


def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f run_data.star")
	parser.add_option("-f",dest="fpar",type="string",metavar="FILE",
		help="run_data.star")

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
        l1 = [x for x in ll1 if len(x.split())>10]
	FORMAT = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d %6d%9.1f%9.1f%8.2f%8.2f%10d%11.4f%8.2f%8.2f\n"
	fout = file("%s_relion2freali_cs_MT.par"%params['fpar'][:-5],"w")
	count = 0


	for i in ll1[:100]:
		if '_rlnImageName' in i:
                        Micro = int(i.split('#')[-1])
                elif '_rlnDefocusU' in i:
                        DefU = int(i.split('#')[-1])
                elif '_rlnDefocusV' in i:
                        DefV = int(i.split('#')[-1])
                elif '_rlnDefocusAngle' in i:
                        DefA = int(i.split('#')[-1])
		elif '_rlnAngleRot' in i:
			Rot = int(i.split('#')[-1])
                elif '_rlnAngleTilt' in i:
                        Tilt = int(i.split('#')[-1])
                elif '_rlnAnglePsi' in i:
                        Psi = int(i.split('#')[-1])
                elif '_rlnOriginXAngst' in i:
                        OriginX = int(i.split('#')[-1])
                elif '_rlnOriginYAngst' in i:
                        OriginY = int(i.split('#')[-1])
                elif '_rlnHelicalTubeID' in i:
                        HexID = int(i.split('#')[-1])

	tt = l1[0].split()
	
	#apix = float(tt[Det-1])/float(tt[Mag-1])*10000
	#apix = 2.68
	#print "apix = %.3f"%apix
        occ = 100
        logP = -6500
        sigma = 0.6
        score = 15.0
        change = 0.0
	mag = 10000.
	MTID_prev = -1
	MT = 0
	frame_prev = -1
	psi_prev = float(tt[Psi-1])
	#print psi_prev

	for i in l1:
		count += 1
                t1 = i.split()
		frame = t1[Micro-1].split('/')[-1]
                df1 = float(t1[DefU-1])
                df2 = float(t1[DefV-1])
                angast = float(t1[DefA-1])
		rot_r = float(t1[Rot-1])
		tilt_r = float(t1[Tilt-1])
                psi_r = float(t1[Psi-1])
		shx = float(t1[OriginX-1])
                shy = float(t1[OriginY-1])
		psi = psi_r
		theta = tilt_r
		phi = rot_r
		MTID = int(t1[HexID-1])
		#if count < 10:
		#	print psi,psi_prev
                #if frame != frame_prev or ((abs(psi-psi_prev) > 6.0 and abs(abs(psi-psi_prev)-180) > 6.0)):
		if MTID != MTID_prev:
                        MT += 1
		MTID_prev = MTID
                frame_prev = frame
		psi_prev = psi
		fout.write(FORMAT%(count,psi,theta,phi,-1*shx,-1*shy,mag,MT,df1,df2,angast,occ,logP,sigma,score,change))
	

	f1.close()
	fout.close()

if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
