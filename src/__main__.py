#!/usr/local/bin/python3.5

import globalvariables as o
import basis
from memtime import memtime
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
# -------------------------------------------------------


# -------------------------------------------------------
# Now do what have to be done!
# -------------------------------------------------------


# -------------------------------------------------------
# read input file param.py and set various global variables
import readinput;
memtime('readinput');# print memory and time info
# start some local variables for use below
nlist = o.nlist;
n, m, mx, Np = o.n, o.m, o.mx, o.Np;
corrtd, matelem = o.corrtd, o.matelem;
groundstate, justenergy = o.groundstate, o.justenergy

# -------------------------------------------------------
niter=0;
for N in nlist:
	o.n = N;
	print('');
	print('========= n = '+str(N)+' ======== '+str(niter+1)+'/'+str(len(nlist))+' ======== ');
	print('');
	# calc symmetrised phonon basis related quantities
	# sets global: Nv1l,Nv2l,Norm1l,Norm2l,Norm3l,map21,map32
	# also sets: n1fsym,n1,n2[,n3],ntot, listn2 [,listn3]
	basis.fbasis(o.n,m,mx,Np);
	memtime('basis');# print memory and time info
	# -------------------------------------------------------
	# calculate Hamiltonian:
	# set global variables: Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
	if niter==0: import hamiltonian; 
	else: reload(hamiltonian)
# import here to use the right values for global variables
	hamiltonian.hamilt();
	memtime('hamiltonian');# print memory and time info
	# -------------------------------------------------------
	if corrtd:
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
			print(" input nstates < 0 or > ntot! so setting nstates = ntot")
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
		if not justenergy:
			if niter==0: import densitymatrices
			else: reload(densitymatrices)
			densitymatrices.cdms();
			memtime('cdms');# print memory and time info

	niter +=1;
memtime('everything');# print memory and time info
print(" ===> Everything completed...")
print(' ')

