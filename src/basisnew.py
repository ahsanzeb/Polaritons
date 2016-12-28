
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
from multiprocessing import Pool
#from pathos.multiprocessing import ProcessingPool as Pool
import decimal
import mapping
import time


#########################################################################
def fbasis(n,m,mx,Np):
	# sets global: Nv1l,Norm1l,Norm2l,Norm3l,map21
	print(' ====> calculating basis states ... ');
	# -------------------------
	if n == 1:
		fbasisn1(n,mx);
		return
	res=[];
	# start lists for one site:
	for i in range(0,m+1):
		res.append([i]);
	# -------------------------
	# use pool or single process? & if pool, best nrest?
	Nbc=1000;
	if n<4  or int(scipy.special.binom(m+n, n)) <= Nbc:
		askpool = 0; nrest = n-1; # all rest iterations in mkfsgnnf()
		print(' a single process .... ')
	else:
		print(' pool of process .... ')
		askpool = 1;
		# ---------------------------------
		# start pool ASAP but to divide the load evenly,
		# start with a large enough set, len(set) ~ 1000*Np
		# increase size of elem of lists to (n-nrest) sites:
		# nrest >= 1 
		res2=[]; nrest = n-1;
		while (len(res) < 1000*Np and nrest > 1):
			for j in res:
				# prepend i in [0,M] for every j with j[0] <= i
				for i in range(j[0],m+1):
					x=[i]+j
					res2.append(x)
			res=res2; res2=[];
			nrest -= 1;
		# print('nrest = ',nrest)
	# ---------------------------------
	# set global array to be used in mkfsgnnf()
	o.res = res;
	o.nrest = nrest;
	# print('after n-nrest: len(res) = ',len(res))
	# -----------------------------
	# why a single function (mkfsgnnf) to do everything:
	# to avoid memory multiplication for pool,
	# the Norms and Nvlist etc are also calculated 
	# with the sets in the same function.
	# (if we first calc sets and then create another pool for the Norms etc calc, the memory usage will increase to Np times, which is bad here because sets (res) contain N elements each so it's not like creating a pool with at diagonalisation time where large matrix H is sparse with far fewer elements in each row than Np.)
	#--------------------------		
	Factlist=[]; Nvlist=[]; Norm1=[]; Norm2=[]; Norm3=[];
	if not askpool:
		listres = [0,len(res)];
		if o.corrtd and o.ld>0:
			Factlist,Nvlist,Norm1,Norm2,Norm3 = mkfsgnnf(listres)
		else:
			Nvlist,Norm1,Norm2,Norm3 = mkfsgnnf(listres)
	elif askpool:
		# listres = mkchunks(Np,len(res));
		# print(' listres = ',listres)
		listres = mkevenloadlist(res,nrest,n,m,Np);
		# print(' even listres = ',listres)
		# print(listres)
		pool=Pool(Np);
		results = pool.map(mkfsgnnf, listres);
		pool.close();
		# combine the arrays from pool
		if o.corrtd and o.ld>0:
			#---------------------------------
			# python list addition is much faster than creating numpy array
			#---------------------------------
			# t0 = time.time();
			# dt=np.dtype('uint32');dt2=np.dtype('float64');
			# sz = int(scipy.special.binom(m+n, n));
			# Factlistnp = np.zeros(sz,	dtype=dt2);
			# Nvlistnp = np.zeros(sz,	dtype=dt);
			# Norm1np = np.zeros(sz,	dtype=dt);
			# Norm2np = np.zeros(sz,	dtype=dt);
			# Norm3np = np.zeros(sz,	dtype=dt);
			# l = 0;
			# for Fctl, Nvl, Nrm1, Nrm2, Nrm3 in results:
			# 	lt = len(Nvl);
			# 	#Factlistnp[l:l+lt]= [float(x) for x in Fctl]
			# 	Nvlistnp[l:l+lt]= Nvl
			# 	Norm1np[l:l+lt]= Nrm1
			# 	Norm2np[l:l+lt]= Nrm2
			# 	Norm3np[l:l+lt]= Nrm3
			# 	l += lt;
			# t1 = time.time();
			#---------------------------------
			for Fctl, Nvl, Nrm1, Nrm2, Nrm3 in results:
				Factlist += Fctl; 
				Nvlist += Nvl;
				Norm1 += Nrm1; Norm2 += Nrm2; Norm3 += Nrm3;	
			# t2 = time.time();
			# print('time for np arrays + fctl float conversion  = ',t1-t0)
			# print('time for python lists  = ',t2-t1)
		else:
			il = 0
			for Nvl, Nrm1, Nrm2, Nrm3 in results:
				Nvlist += Nvl;
				Norm1 += Nrm1; Norm2 += Nrm2; Norm3 += Nrm3;
				il += len(Nrm1);
				print('len(Nrm1) = ',len(Nrm1))
			print('tot = ',il)
		del results; 
	#--------------------------
	n2 = int(scipy.special.binom(m+n-1, n-1));# basis size for N-1 sites
	Norm2 = Norm2[0:n2];
	if n>2:
		n3 = int(scipy.special.binom(m+n-2, n-2)); # basis size for N-2 sites
		Norm3 = Norm3[0:n3]
	else:
		n3 = 0; Norm3 = [];
	#--------------------------
	# set global arrays:
	o.Fact1l = Factlist;
	o.Nv1l = Nvlist;
	o.Norm1l, o.Norm2l, o.Norm3l = Norm1,Norm2,Norm3
	#--------------------------
	n1fsym = int(scipy.special.binom(m+n, n)); # old thing with name n1
	n1 = n1fsym + (mx-m)*n2; # extended basis size, photon sector
	# total number of basis states
	ntot = n1 + (mx+1)*n2;	
	o.n1fsym,o.n1,o.n2,o.n3,o.ntot = n1fsym,n1,n2,n3,ntot; 
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
	#print(' calculating mapping ... ')
	#mapping.getmapp(n,m,Np);# calculate and set o.map21
	return
#-----------------------------------------------
# make full sets and get nv,norms,fac
def mkfsgnnf(indlist):
		# uses res from parent function
		# repeat nrest times to get sets for n sites.
		nrest = o.nrest;
		res = o.res;
		n=o.n; m = o.m;
		#print('res = ',res);
		#print('indlist = ',indlist);
		res2=[]; 
		i1 = indlist[0]; i2 = indlist[1];
		resx = res[i1:i2]
		for dummy in range(nrest):
			for j in resx:
				# prepend i in [0,M] for every j with j[0] <= i
				for i in range(j[0],m+1):
					x=[i]+j
					res2.append(x)
			resx=res2
			res2=[];
		# get no of phonons list, Normalisation list, etc.
		if o.corrtd and o.ld>0:
			# Fctl,Nvl,Nrm1,Nrm2,Nrm3 = getNvNormFacl(resx,n,m);
			return getNvNormFacl(resx,n);# Fctl, Nvl, Nrm1, Nrm2, Nrm3
		else:
			# Nvl,Nrm1,Nrm2,Nrm3 = getNvNorml(resx,n,m);			
			return getNvNorml(resx,n);#Fctl, Nvl, Nrm1, Nrm2, Nrm3
		return
#-----------------------------------------------
# also calc factorials for corrtd psi0 with ld>0
def getNvNormFacl(resx,n):
		Nvlist =[]; Factlist=[];
		Norm1list=[]; Norm2list=[]; Norm3list=[];
		m0=math.factorial(n) # normalisation factor
		cont=decimal.Context(prec=15, Emax=999, clamp=1)
		# for n=2, Norm3list is not meant to be used!
		for j in resx:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			# if this set index lies in range 0-n2,
			# the first elem is always 0. 
			# to use this function alone, we comput of 
			# m2 and m3 for all sets in res, and after 
			# combining all the results from pool,
			# we take only valid data for Norm2l and NOrm3l,
			# i.e., only from fisrt n2, n3 sets, and discard the rest.	
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
		return Factlist, Nvlist, Norm1list, Norm2list, Norm3list
# --------------------------------------------
def getNvNorml(resx,n):
		Nvlist =[];
		Norm1list=[]; Norm2list=[]; Norm3list=[];
		m0=math.factorial(n) # normalisation factor
		# for n=2, Norm3list is not meant to be used!
		for j in resx:
			nv=np.sum(j) # total no of phonons
			Nvlist.append(nv)
			k=stats.itemfreq(j)
			# if this set index lies in range 0-n2,
			# the first elem is always 0. 
			# to use this function alone, we comput of 
			# m2 and m3 for all sets in res, and after 
			# combining all the results from pool,
			# we take only valid data for Norm2l and Norm3l,
			# i.e., only from fisrt n2, n3 sets, and discard the rest.	
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
		return Nvlist, Norm1list, Norm2list, Norm3list
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
# create list of indices ranges that lead roughly to even load distribution in pool map.
def mkevenloadlist(res,nrest,n,m,Np):
	eloadl = [];
	nsets =np.zeros(len(res));
	l = 0;
	for j in res:
		j0 =j[0]; mm = m-j[0];
		nn = int(scipy.special.binom(mm+nrest,nrest));
		nsets[l]=nn; l += 1;
	nntot = np.sum(nsets);
	dn = nntot/Np;
	#print('nsets = ',nsets)
	#print('dn = ',dn)
	i = 0; i1 = 0; sm = 0; v2 = 1*dn; ichunk = 0;
	while i < l-1 and ichunk < Np-1:
		sm += nsets[i];
		if sm >= v2:
			if i == i1:
				i2 = i+1;
				#print(' i=i1 ')
			elif abs(sm-v2)<abs(sm-nsets[i]-v2):
				i2 = i;
				#print(' < ')
			else:
				i2 = i-1;
				#print(' > ')				
			eloadl.append([i1,i2]); ichunk += 1;
			i1 = i2; v2 += dn;	
			#i = i2 +1;
		i += 1;	
	#print(' last one =',[i1,l])	
	eloadl.append([i1,l]); # last chunk
	return eloadl
#-----------------------------------------------


