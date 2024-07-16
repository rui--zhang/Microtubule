#!/usr/bin/env python
import math
import os

f1 = file('runparfile',"w")
#f2 = file('f-list-par.txt',"w")
#f3 = file('f-list-shft.txt',"w")
#curpath = os.getcwd().replace('/jetstor/CHOME/cluster/ruiz','/home/ruiz/cluster')

f2 = file('mparameters')
l2 = f2.readlines()
iter = int(l2[-4].split()[-1])
print "iter = %d"%iter

start = int(l2[-3].split()[-1])
last = int(l2[-2].split()[-1])
print "first ptcl = %d"%start
print "last ptcl = %d"%last

#nodes = 1
# on 16G nodes, each takes 13.8% memory, 13.8% x 6 = 82.8%
nproc = 50
ncls = 1
print "ncls = %d"%ncls

incr = math.ceil((last-start)/nproc)+1   # increment
startn = start - incr
lastn = incr

while(lastn <= last):
	startn = startn + incr 
	lastn = startn + incr - 1
	x = min(lastn,last)
	f1.write("./mult_hrefine_n.com %d %d %d %d&\n"%(startn,x,iter,ncls))
	f1.write("sleep 1\n")

f1.close()
f2.close()

os.system('head -70 runparfile >runparfile1')
os.system('tail -70 runparfile >runparfile2')
os.system('chmod a+x runparfile*')
