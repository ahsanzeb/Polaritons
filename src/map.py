
# construct the map of indices from N-1 site to N-site.

import scipy
import numpy as np

# starting with x0,y0, the recursive function
# makes the indices in upper & lower triangular, 
# and the left down column, and 
# gives coordinates (& starting value)for next block top left, .... repeats, until the bottom right corner is reached.

#-----------------------------
def mkarglist(list1):
	list2l=[]
	for x in list1:
		list2l+=list(range(1,x+1))[::-1];
	return list2l
#-----------------------------
def mapn2(m):
	k = 0;
	dt = np.dtype('uint32');
	Map = np.zeros((m+1,m+1),dtype=dt);
	for i in range(m+1):
		for j in range(i,m+1):	
			Map[i,j] = int(k);
			Map[j,i] = int(k);
			k += 1;
	return Map
#-----------------------------
def triangles(colist,iin,k,n,m):
		nrows = (iin+1)*iin//2;
		dt = np.dtype('uint32');
		Map = np.zeros((nrows,m+1),dtype=dt);
		# -------------------------
		# fill all left side columns for this iin
		i = 0;
		for x in colist:
			Map[:,i] = range(x,x+nrows);
			i += 1;
		# -------------------------		
		yi = 0; i=0; 
		for i in range(iin):
			i0 = i+m-iin+1;
			yf = (i+1)*iin - i*(i+1)//2;
			for jj in range(yi,yf):
				j = jj - yi;
				for mjj in range(i0+j,m+1):
					mj = mjj -i0;
					Map[jj,mjj] = int(k);
					Map[yi+mj,i0+j] = int(k);
					if j ==0 and mjj==m: km = k;
					k += 1;
			Map[yf:,i0] = range(km+1,k)
			yi = yf;
			i += 1;
		return k, Map
#-----------------------------
def getmap(n,m):
	mapfull = mapn2(m);	# n = 2 case
	karg = mapfull[-1,-1] +1;
	arglist=[m]; # starting argument list
	for nn in range(3,n+1):
		for iarg in arglist:
			colist = mapfull[-1,0:m-iarg+1]+1;
			karg,mapl = triangles(colist,iarg,karg,nn,m);
			mapfull=np.concatenate((mapfull,mapl));
		arglist=mkarglist(arglist); # make argument list
	return mapfull
#-----------------------------

if 1:
	n =10; m=5;
	mapfull = getmap(n,m)


