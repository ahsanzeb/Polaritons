#!/usr/local/bin/python3.5

import globalvariables as o
import basis
from memtime import memtime
from auxfunc import prntmsg,createoutdir,setNp


import time, sys
from importlib import reload
import importlib

# -------------------------------------------------------
# set the input file name for readinput module:
if len(sys.argv) > 1: o.inpfile = sys.argv[1];
else: o.inpfile = 'param.py';
# -------------------------------------------------------
# welcome message 
print(" ")
print(" *** *** *** Welcome to Polaritons *** *** *** ")
print(" Current date & time " + time.strftime("%c"))
import bag
# -------------------------------------------------------
# read input file param.py and set various global variables
import readinput;
memtime('readinput');# print memory and time info
# start some local variables for use below
nlist = o.nlist; 
mlist = o.mlist;
n, m, mx, Np0 = o.n, o.m, o.mx,o.Np;
corrtd, matelem = o.corrtd, o.matelem;
zeroTPL = o.zeroTPL;
groundstate, justenergy = o.groundstate, o.justenergy
photonfraction = o.photonfraction;

diffoutdir = o.diffoutdir;
if zeroTPL:
	print(" PL: undisplaced basis will be used... ")
	o.ld = 0;
old = o.ld; # to reset ld value for n!=1
lamlist = o.lamlist;
uselamlist = o.uselamlist;
wclist=o.wclist;
wvlist=o.wvlist;
# -------------------------------------------------------
niter=0; lnlist = len(nlist);
for nn in nlist:
	n = int(nn);
	o.n = n;
	m, mx = mlist[niter];
	o.m, o.mx = m, mx;

	print('lambin0 = ', o.lambin0)
	
	if uselamlist:
		o.lambin0 = [lamlist[niter]]; # to cheat
		print(' lam = ',o.lambin0)

	#o.lambin0 = [0.01,0.02,0.03,0.04];
	#print('lambin0 new = ', o.lambin0)

	if len(wclist)==len(nlist):
		o.wc = wclist[niter]; print(' wc = ',o.wc)
	if len(wvlist)==len(nlist):
		o.wv = wvlist[niter]; print(' wv = ',o.wv)

	# -------------------------
	# print iter message
	prntmsg(n,m,mx,niter,lnlist); 
	# create output directory
	o.dumy = createoutdir(n,diffoutdir);
	if n==1 or o.nlarge: #or m==0:
		o.ld = 0; # undisplaced basis are used
	else:
		o.ld = old;
	
	#o.ld = 0.1*niter;
	print('o.ld = ',o.ld)

	print('main: o.nlarge =', o.nlarge)
	# set Np <= n3
	Np = setNp(n,m,Np0); o.Np=Np;
	# -------------------------
	# calc symmetrised phonon basis related quantities
	# sets global: Nv1l,Nv2l,Norm1l,Norm2l,Norm3l,map21,map32
	# also sets: n1fsym,n1,n2[,n3],ntot, listn2 [,listn3]
	if(o.nlarge):
		basis.fbasisNormsRatio(o.ndummy,m,mx,Np);
	else:
		basis.fbasis(n,m,mx,Np);

	memtime('basis');# print memory and time infox
	# -------------------------------------------------------
	# calculate Hamiltonian:
	# set global variables: Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
	if niter==0: import hamiltonian; 
	else: reload(hamiltonian)
	# import here to use the right values for global variables
	hamiltonian.hamilt(o.nlarge);
	memtime('hamiltonian');# print memory and time info	
	# -------------------------------------------------------
	if corrtd:
		if zeroTPL:
			if niter==0: import plmatelem; 
			else: reload(plmatelem)
			plmatelem.plmatelem();
		else:
			# absorption: using time evolution:
			# calculates and saves: corr,ft, analyt. G
			if niter==0: import correlation; 
			else: reload(correlation)
			correlation.corrft();
		memtime('corrft');# print memory and time info
	# -------------------------------------------------------	
	if matelem:
		# absorption: matrix elements:
		if o.nstates < 0 or o.nstates > o.ntot:
			o.nstates = o.ntot;
			o.nstates = o.n1	;
			print(" input nstates < 0 or > ntot! so setting nstates = ntot: nstates = ",o.nstates)
		if niter==0: import matrixelements;
		else: reload(matrixelements)
		matrixelements.fmatelem();
		memtime('matelem');# print memory and time info
	# -------------------------------------------------------
	if groundstate:
		# ground state calculations:
		if niter==0: import gssolver;
		else: reload(gssolver)
		gssolver.gssolve();
		memtime('gssolve');# print memory and time info
		# write energy output file
		if niter==0: import energy;
		else: reload(energy)
		if photonfraction:
			if niter==0: import photon
			else: reload(photon)
			photon.photfrac()
		if not justenergy:
			if niter==0: import densitymatrices
			else: reload(densitymatrices)
			densitymatrices.cdms();
			memtime('cdms');# print memory and time info

	niter +=1;
memtime('everything');# print memory and time info
print(" ===> Everything completed...")
print(' ')

