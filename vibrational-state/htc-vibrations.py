from qutip import *
from numpy import linspace, meshgrid, sqrt
from scipy.sparse.linalg import lobpcg
from scipy.sparse import identity
import numpy as np

import os,sys
os.makedirs('fullHilbert/data/n-all', exist_ok=True)

# path to input parameter file
path_param = os.getcwd();
# add this path to python search path
sys.path.insert(0, path_param)




wv=0.2;
wr=1;
lamlist=[0.5,1,1.5,2,2.5];

N=5;
m=N;
M=5;


#---------------------------------------
def mkHb(N,M):
	sp=	create(2); sm=	destroy(2);
	I2=eye(2); spm=sp*sm;
	IM=eye(M); Xb=create(M)+destroy(M);
	Hb= ;# initialise Hb as qobj()
	for i in range(N):
		x=[];
		# electronic 2-LS
		for j in range(N):
			if j == i:
				x.append(spm);
			else:
				x.append(I2);
		# vibrational M-LS
		for j in range(N):
			if j == i:
				x.append(Xb);
			else:
				x.append(IM);
		# vibronic part of H
		Hb += tensor(x);
	return Hb
#---------------------------------------
def getNup(N,m):
	I2=qeye(2); spm=create(2)*destroy(2);
	Hb= ;# initialise Hb as qobj()
	for i in range(N):
		x=[];
		for j in range(N):
			if j == i:
				x.append(spm);
			else:
				x.append(I2);
		Hb += tensor(x);
	return m-Hb;
#---------------------------------------
def mkHR(N,m):
	sp=	create(2); I2=eye(2);
	Nup = m - mkN(N,2);
	Nup = Nup.sqrtm();
	Hb= ;# initialise Hb as qobj()
	for i in range(N):
		x=[];
		# electronic 2-LS
		for j in range(N):
			if j == i:
				x.append(sp);
			else:
				x.append(I2);
		Hb += tensor(x);
	# muplitply the sqrt(Nphoton) factors column-wise
	Hb = Hb*Nup		
	return Hb
#---------------------------------------
def mkI(N,M):
	Ix=eye(M);
	X=[];
	for j in range(N):
		X.append(Ix);
	return tensor(X);
#---------------------------------------
def mkN(N,M):
	I2=qeye(M); spm=create(M)*destroy(M);
	Hb= ;# initialise Hb as qobj()
	for i in range(N):
		x=[];
		for j in range(N):
			if j == i:
				x.append(spm);
			else:
				x.append(I2);
		Hb += tensor(x);
	return Hb;
#---------------------------------------

if m<N:
	print("m<N not allowed! stop!");
	exit();
	
Hr = mkHR(N,m) * mkI(N,M);
Ielec = mkI(N,2);
Hb = Ielec * mkHb(N,M);
Hv = Ielec * mkN(N,M);

H0 = wv*Hv + wr*Hr;
for lam in lamlist:
 H = H0 + wv*lam*Hb;
 psi= H.groundstate();
 



