

# calculate Nv{1,2,3}List, mapping21 and mapping32
# Norm1list,Norm2list,Norm3list,Hgcoor21,Hgcoor32 

#########################################################################
# Ahsan Zeb, 03 Jun 2016
#########################################################################
from scipy import stats
import numpy as np
import math
import sys

import globalvariables as o

#########################################################################
# IN: n = no. of sites, m = max no. of phonon per site
# OUT: Norm1list,Norm2list,Norm3list,Hgcoor21,Hgcoor32
def fbasis(n,m,mx,Np):
	# sets global: Nv1l,Nv2l,Norm1l,Norm2l,Norm3l,map21,map32
	print(' ====> calculating basis states ... ')
	res=[];
# -------------------------
# start lists for one site:
	for i in range(0,m+1):
		res.append([i]);
# -------------------------
# increase size of elem of lists to n-2 sites:
#	dnt execute for n<4 cases: range(2,n-1)
	res2=[];
	for dummy in range(2,n-1): 
		for j in res:
			# prepend i in [0,M] for every j with j[0] <= i
			for i in range(j[0],m+1):
				x=[i]+j
				res2.append(x)
		res=res2
		res2=[];

# Normalisation list
# only if n > 2
	if ( n > 2):
		Norm3list = getNorm(res,n-2)
		ntot3 = len(Norm3list)
	else:
		Norm3list=None; ntot3=None

# ---------------------------------------
# add a site, to get basis for N-1 site and
# calculate links32, Hgcoor32 for N-2 to N-1 sites
# these will be used to calculate full mapping matrix
# for N-2 site to N-1 site basis
# ---------------------------------------
# only if n > 2
	if ( n > 2):
		res,links32,Hgcoor32 = addsite(res,m);
		# ---------------------------------------
		# calculate missing coordinate mapping for 2to1:
		Hgcoor32mis = fHgcoor2(links32,ntot3,m)
		del links32
		# ---------------------------------------
		# calculate full mapping matrix for 32:
		Hgcoor32 = getfullmap(Hgcoor32,Hgcoor32mis,ntot3,m)
		del Hgcoor32mis
	else:
		Hgcoor32=None
		
#	print(Hgcoor32)
# ---------------------------------------
# No of phonons list, Normalisation list
	Nv2list,Norm2list = getNvNorm(res,n-1)
	ntot2 = len(Norm2list)
	
# ---------------------------------------
# add a site, to get basis for N site and
# calculate links21, Hgcoor21 for N-1 to N sites
	res,links21,Hgcoor21 = addsite(res,m);
# calculate missing coordinate mapping for 2to1:
	Hgcoor21mis = fHgcoor2(links21,ntot2,m)
	del links21
# calculate full mapping matrix for 21:
	Hgcoor21 = getfullmap(Hgcoor21,Hgcoor21mis,ntot2,m)
	del Hgcoor21mis
#	print(Hgcoor21)
# No of phonons list, Normalisation list
	Nv1list,Norm1list = getNvNorm(res,n)
	ntot1 = len(Norm1list)

# 


# ---------------------------------------
#	print(" n1, n2, n3 = ", ntot1,ntot2,ntot3)

	o.Nv1l, o.Nv2l = Nv1list, Nv2list
	o.Norm1l, o.Norm2l, o.Norm3l = Norm1list, Norm2list, Norm3list
	o.map21,o.map32 = Hgcoor21, Hgcoor32

	# -------------------------------------------------------
	# also calculate and set: n1fsym,n1,n2[,n3],ntot, listn2 [,listn3]
	# -------------------------------------------------------
	#from functions import countbasis, mkpoollists
	# get counts on various basis sectors
	# and set global variables: n1fsym,n1,n2[,n3],ntot
	#countbasis();
	n1fsym=len(Norm1list); # old thing with name n1
	n2=len(Norm2list); # No. of symmetrised basis for |0c,1x>
	# new n1: so that we dont have change shifts by n1 everywhere
	n1 = n1fsym + (mx-m)*n2;
	# total number of basis states
	ntot = n1 + (mx+1)*n2;
	if (n>2):
		n3=len(Norm3list); # basis size for N-2 sites
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
	return

#########################################################################

# calculate mapping
	# map21 = []; map32 = []; 


####################OLD COMMENTS###############
#   recording linksjj2i {njj2i,pntri} list.
#  njj2i = Number of basis in sets (basi) formed from a given basis in sets2 (basjj); 
#  pntri = pointer to start index of corresponding basi.
####################OLD COMMENTS###################
# calculate matrix elements of Hg missed in Hgcoor
# using a more direct method
# nj2i is a list/nparray of [njj2i,pnti]. 
# pool map on mj, bcs it's simple if we have full loop over basjj.
########################OLD COMMENTS############
# in notes: m stands for m0, m0=m+1 ==> z = m+2-mj
# Here m is m, and extra +1 is dropped bcs mj starts frm 0
# so z = m+1-mj; 
#########################################################################
# nj2i is list of links between the basis sets with sites info in 'inf'
def	fHgcoor2(nj2i,n22,m):
	Hgcmis = [];
	for mj in range(0,m+1):
		z = m+1-mj;
		jj = 0; ii = 0;
		# jj = 0 to n2-1; n2 items
		while (jj < n22): 
			# find start of a contributing patch
			while(jj < n22 and nj2i[jj][0] >= z): 
				jj += 1;
			#-------end of non-contrib patch-----------
			# down-left-next DLN rule:
			ii = nj2i[jj-1][1] + 1 ; # corresponding basi index;
			# calculate hg's mat elem in this patch
			while(jj < n22 and nj2i[jj][0] < z):
				Hgcmis.append([jj,mj,ii]) # note the order
				ii += 1; # add 1
				jj += 1; # add 1
		#--------end of contrib patch----------
	#----- considered complete basjj set ----
	return Hgcmis
#########################################################################
# ---------------------------------------
# No of phonons list, Normalisation list
# ---------------------------------------
def getNvNorm(resx,nsit):
	Nvlist =[]; Normlist=[];
	m0=math.factorial(nsit) # normalisation factor
	for j in resx:
		nv=np.sum(j) # total no of phonons
		Nvlist.append(nv)
		k=stats.itemfreq(j)
		m1=m0
		for kk in k:
			# normalisation factor
			rs=math.factorial(kk[1]) 
			m1=m1/rs # divide by second entry, i.e., frequency
		Normlist.append(m1); #append(int(m1))
	return Nvlist, Normlist
#########################################################################
# ---------------------------------------
# Normalisation list
# ---------------------------------------
def getNorm(resx,nsit):
	Normlist = [];
	m0=math.factorial(nsit) # normalisation factor
	for j in resx:
		k=stats.itemfreq(j)
		m1 = m0 
		for kk in k:
			rs=math.factorial(kk[1]) 
			m1=m1/rs # divide by second entry, i.e., frequency
		Normlist.append(m1);
	return Normlist
#########################################################################
# add a site, keep record of links of old to new basis set in hgcor & links
# ---------------------------------------
# make sets with N elements
# keep track of indeces of prev sets, new sets, and the elem added
# This is to calculate most* of the matrix elements of Hg without finding
# these quantities again, thus saving time.
# * hgcoor corresponding to sets; mj>=sets2[0]
def addsite(resx,m):
	res2 = [];
	links=[]; Hgcor=[]; n21=0; q0=0;
	p=-1; q=-1;
	# take every set
	for j in resx: 
		p=p+1; #  index of basis for N-1 site system
		# prepend i in [0,M] for every j with j[0] <= i
		for i in range(j[0],m+1):
			q=q+1; #  index of basis for N site system
			x=[i]+j
			# make updated sets
			res2.append(x)
			# coor of Hg non-zero mat elem
			Hgcor.append([p,i,q])
			n21 += 1;
		links.append([n21,q0]);
		n21 = 0; q0=q+1;
	return res2, links, Hgcor
#########################################################################
# calculate full mapping matrix:
def getfullmap(Hgcor,Hgcormis,nt,m):
	# join lists
	Hgcor += Hgcormis
 	# Sort by first elements of the inner lists, basjj
	Hgcor.sort()
 	# make a list of only last elements, basi
	Hgcor=np.array(Hgcor) # convert to np.array
 	# Hgcor is ordered wrt first and second col.
 	# just keep third col. & reshape to make a matrix
	Hgcor = Hgcor[:,2].reshape((nt,m+1))
	return Hgcor
#########################################################################

# write sets in frequency format
def freq(aset):
	asf=np.zeros(m+1);
	for elem in aset:
		asf[elem] += 1
	return asf
#########################################################################






import globalvariables as o


#****************************************************************
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
#****************************************************************

