
import globalvariables as o

# calculate Nv{1,2,3}List, mapping21 and mapping32
# Norm1list,Norm2list,Norm3list,Hgcoor21,Hgcoor32 

#########################################################################
# Ahsan Zeb, 03 Jun 2016
#########################################################################
import scipy
from scipy import stats
import numpy as np
import math
import sys

import decimal
import mapping

#########################################################################
def fbasis(n,m,mx,Np):
	# sets global: Nv1l,Norm1l,Norm2l,Norm3l,map21
	print(' ====> calculating basis states ... ');
	# -------------------------
	if n == 1:
		fbasisn1(n,mx);
		return
	# -------------------------
	res=[];
	# start lists for one site:
	for i in range(0,m+1):
		res.append([i]);
	# -------------------------
	# increase size of elem of lists to n sites:
	res2=[];
	for dummy in range(1,n): 
		for j in res:
			# prepend i in [0,M] for every j with j[0] <= i
			for i in range(j[0],m+1):
				x=[i]+j
				res2.append(x)
		res=res2
		res2=[];

# No of phonons list, Normalisation list
	if o.corrtd and o.ld>0:
		getNvNormFac(res,n,m); # sets global arrays
	else:
		getNvNorm(res,n,m); # sets global arrays

	n1fsym = int(scipy.special.binom(m+n, n)); # old thing with name n1
	# No. of symmetrised basis for |0c,1x>
	n2 = int(scipy.special.binom(m+n-1, n-1)); 
	# new n1: so that we dont have change shifts by n1 everywhere
	n1 = n1fsym + (mx-m)*n2;
	# total number of basis states
	ntot = n1 + (mx+1)*n2;
	
	if (n>2):
		n3 = int(scipy.special.binom(m+n-2, n-2)); # basis size for N-2 sites
		o.n1fsym,o.n1,o.n2,o.n3,o.ntot = n1fsym,n1,n2,n3,ntot; 
	o.n1fsym,o.n1,o.n2,o.ntot = n1fsym,n1,n2,ntot;
	# -------------------------------------------------------
	# get chunk lists for multiprocessing pool map
	# sets global variables: listn2 [,listn3]
	# mkpoollists(); 
	# -------------------------------------------------------
	o.listn2 = mkchunks(Np,n2);
	if (n>2):
		o.listn3 = mkchunks(Np,n3);
	# -------------------------------------------------------
	print(" n_photon, n_exciton, ntot = ", n1, n2*(mx+1), ntot)
	print(' calculating mapping ... ')
	mapping.getmap(n,m);# calculate and set o.map21
	return

#-----------------------------------------------
# seperate copies of code operating in different ranges
# instead of if statements to seperate these cases.
# hopefully it would be more efficient this way.
#-----------------------------------------------
# also calc factorials for corrtd psi0 with ld>0
def getNvNormFac(res,n,m):
	Nvlist =[]; Factlist=[];
	Norm1list=[]; Norm2list=[]; Norm3list=[];
	m0=math.factorial(n) # normalisation factor
	cont=decimal.Context(prec=15, Emax=999, clamp=1)
	n3 = 0;
	n2 = int(scipy.special.binom(m+n-1, n-1));
	# obtain everything for n-2 site case:
	if n >2:
		n3 = int(scipy.special.binom(m+n-2, n-2));
		for j in res[0:n3]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			rzero = k[0,1]; # number of sites with 0 vibrations
			m1=m0; fc0 = decimal.Decimal.from_float(1);
			for kk in k:
				# n0=factorial of occupation numbers
				n0=math.factorial(kk[0]);
				n0 = decimal.Decimal.from_float(n0);
				fc0 *= cont.power(n0,int(kk[1]));
				# normalisation factor
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Factlist.append(fc0);
			Norm1list.append(m1);
			m2 = m1*rzero/n;
			Norm2list.append(m2);
			m3 = m2*(rzero-1)/(n-1);
			Norm3list.append(m3);
	# obtain everything for n-1 site case:
	for j in res[n3:n2]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			rzero = k[0,1]; # number of sites with 0 vibrations
			m1=m0; fc0 = decimal.Decimal.from_float(1);
			for kk in k:
				# n0=factorial of occupation numbers
				n0=math.factorial(kk[0]);
				n0 = decimal.Decimal.from_float(n0);
				fc0 *= cont.power(n0,int(kk[1]));
				# normalisation factor
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Factlist.append(fc0);
			Norm1list.append(m1);
			m2 = m1*rzero/n;
			Norm2list.append(m2);
	# obtain everything for n site case:
	for j in res[n2:]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			m1=m0; fc0 = decimal.Decimal.from_float(1);
			for kk in k:
				# n0=factorial of occupation numbers
				n0=math.factorial(kk[0]);
				n0 = decimal.Decimal.from_float(n0);
				fc0 *= cont.power(n0,int(kk[1]));
				# normalisation factor
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Factlist.append(fc0);
			Norm1list.append(m1);
	# set global arrays:
	o.Fact1l = Factlist;
	o.Nv1l = Nvlist;
	o.Norm1l, o.Norm2l, o.Norm3l = Norm1list, Norm2list, Norm3list
	return 
# --------------------------------------------
def getNvNorm(res,n,m):
	Nvlist =[];
	Norm1list=[]; Norm2list=[]; Norm3list=[];
	m0=math.factorial(n) # normalisation factor
	n3 = 0;
	n2 = int(scipy.special.binom(m+n-1, n-1));
	# obtain everything for n-2 site case:
	if n >2:
		n3 = int(scipy.special.binom(m+n-2, n-2));
		for j in res[0:n3]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			rzero = k[0,1]; # number of sites with 0 vibrations
			m1=m0;
			for kk in k:
				# normalisation factor
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Norm1list.append(m1);
			m2 = m1*rzero/n;
			Norm2list.append(m2);
			m3 = m2*(rzero-1)/(n-1);
			Norm3list.append(m3);
	# obtain everything for n-1 site case:
	for j in res[n3:n2]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			rzero = k[0,1]; # number of sites with 0 vibrations
			m1=m0;
			for kk in k:
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Norm1list.append(m1);
			m2 = m1*rzero/n;
			Norm2list.append(m2);
	# obtain everything for n site case:
	for j in res[n2:]:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			m1=m0;
			for kk in k:
				rs=math.factorial(kk[1]) 
				m1=m1/rs # divide by second entry, i.e., frequency
			Norm1list.append(m1);
	# set global arrays:
	o.Nv1l = Nvlist;
	o.Norm1l, o.Norm2l, o.Norm3l = Norm1list, Norm2list, Norm3list
	return 
#-----------------------------------------------
#  n =1 case:
#-----------------------------------------------
def fbasisn1(n,mx):
	Nv1list=[]; Norm1list=[];
	for i in range(mx+1):
		Nv1list.append(i);
		Norm1list.append(1);
	o.Nv1l = Nv1list;
	o.Norm1l = Norm1list;
	o.n1 = mx+1; o.ntot = 2*(mx+1);
	return
#-----------------------------------------------
# make even chunks
# divide range(0,nnn) in Npp (= # of processes) intervals with boundaries [i1,i2] in array listn
def mkchunks(Npp,nnn):
	n21 = int(nnn/Npp); rem2 = nnn % Npp;
	listn=[]; i1=0;
	for i in range(rem2):
		i2= i1 + n21 +1
		listn.append([i1,i2])
		i1=i2
	for i in range(rem2,Npp):
		i2= i1 + n21
		listn.append([i1,i2])
		i1=i2
	return listn
#-----------------------------------------------

