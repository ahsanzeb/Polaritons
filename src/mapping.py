
# construct the map of indices from N-1 site to N-site.
import globalvariables as o

import scipy
import numpy as np
from multiprocessing import Pool

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
def triangles(colist,iin,k,m):
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
	#-----------------------------	
	mapfull = mapn2(m);	# n = 2 case
	karg = mapfull[-1,-1] +1;
	arglist=[m]; # starting argument list
	for nn in range(3,n+1):
		for iarg in arglist:
			colist = mapfull[-1,0:m-iarg+1]+1;
			karg,mapl = triangles(colist,iarg,karg,m);
			mapfull=np.concatenate((mapfull,mapl));
		arglist=mkarglist(arglist); # make argument list
	o.map21 = mapfull;
	#print(np.shape(mapfull))
	#print(np.prod(np.shape(mapfull)))
	return 
#-----------------------------




#----------------------------------------------------------
# parallel versions:
#----------------------------------------------------------
def trianglesp(p):
	m = o.m;
	#print('o.eloadlmap = ',o.eloadlmap)
	istart, iend = o.eloadlmap[p];
	nsize = o.nblksizemap[p];
	#print('o.eloadlmap[p]=',o.eloadlmap[p])
	#print('o.nblksizemap = ',o.nblksizemap)
	#print('o.strtindsmap = ',o.strtindsmap)
	#print('o.argxmap = ',o.argxmap)
	#print(' ***********************nsize = ',nsize)
	dt = np.dtype('uint32');
	#print(' ----------   m+1 = ',m+1)
	Map = np.zeros((nsize,m+1),dtype=dt);
	yo = 0;
	for ind in range(istart,iend):
		iin = o.argxmap[ind];
		k = o.strtindsmap[ind] + 1;
		colist = o.colsmap[ind]; #cols is a list of np arrays?

		# -------------------------
		# fill all left side columns for this iin
		nrows = (iin+1)*iin//2;
		i = 0;
		for x in colist:
			#print('nrows, x,yo,yo+nrows = ',nrows, x,yo,yo+nrows)
			Map[yo:yo+nrows,i] = range(x+1,x+nrows+1);
			i += 1;
		# -------------------------		
		yi = yo; i=0; 
		for i in range(iin):
			i0 = i+m-iin+1;
			yf = yo + (i+1)*iin - i*(i+1)//2;
			for jj in range(yi,yf):
				j = jj - yi;
				for mjj in range(i0+j,m+1):
					mj = mjj -i0;
					Map[jj,mjj] = int(k);
					Map[yi+mj,i0+j] = int(k);
					if j ==0 and mjj==m: km = k;
					k += 1;
			Map[yf:yf+k-km-1,i0] = range(km+1,k)
			yi = yf;
			i += 1;
		yo += nrows;	
	return Map
#-----------------------------





#---------------------------------------------------------
# getmapargs(n,m) calculates the arguments for trianglesp
#---------------------------------------------------------
def getmapargs(n,m,Np):
	# -----------------------
	strtinds = []; args = []; cols =[]; trisize = [];
	# starting values for n=3 blocks
	i = (m+1)*(m+2)//2 - 1; # end index of n=2 block
	arg = [m];
	col = np.array([m]);
	# -----------------------
	# get things for for n>3 
	#--------------------
	argx = [];
	for nn in range(3,n+1):
		for ms in arg:
			# use previously calc values
			args.append(ms);
			strtinds.append(i);
			col = col[0:m-ms+1];
			cols.append(col);
			lt = ms*(ms+1)//2;
			trisize.append(lt); 
			# create values for next ms iteration
			col = np.append(col+lt,i+lt);
			i += ms*(ms+1)*(ms+2)//6;
			# create list for next n
			lx = list(range(1,ms+1))[::-1];
			argx += lx;				
		arg = argx; argx=[];

	#-----------------------------
	# set global arrays
	o.argxmap = args;
	o.strtindsmap = strtinds;
	o.colsmap = cols;

	#-----------------------------
	# now make evenload list
	trisize = np.array(trisize);
	# if len(trisize) < 10*Np:
	# 	Np = 1;
	# 	o.Np=Np;
	# 	print(' len(trisize) < 10*Np so setting Np = 1')
	eloadl, nblksize = mkevenloadlist(trisize,Np)
	#print(eloadl)
	# set global arrays: to be used in trianglesp()
	o.eloadlmap = eloadl;
	o.nblksizemap = nblksize;
	#print('XXXX eloadl = ',eloadl)
	#print('XXXX nblksize = ',nblksize)
	return 


#---------------------------------------
# getmapp(n,m,Np) calculates the map
#---------------------------------------
def getmapp(n,m,Np):
	#-----------------------------
	# use single process for small systems:
	if int(scipy.special.binom(m+n-1, n-1))<10000:
		getmap(n,m);
		return
	#--------------------------------
	# calculate the arguments for the 
	# recursive function trianglesp
	# and evenload list for pool.map
	getmapargs(n,m,Np);
	# Np = o.Np; # getmapargs can reset Np
	#--------------------------------
	# now start pool and send jobs
	#--------------------------------
	pool=Pool(Np);
	results = pool.map(trianglesp, range(Np));
	pool.close();
	# join results:
	dt=np.dtype('uint32');
	nsize = int(scipy.special.binom(m+n-1, n-1));
	#print('nsize = ',nsize)
	mapfull = np.zeros((nsize,m+1), dtype=dt);
	mapl = mapn2(m);	# n = 2 case
	nrows = mapl.shape[0];
	mapfull[0:nrows,:] = mapl;
	y1 = nrows;
	for mapl in results:
		lt = mapl.shape[0];
		#print('joining: mapl.shape = ',mapl.shape)
		#print('joining: mapfull.shape = ',mapfull.shape)
		y2 = y1 + lt;
		#print('y1,y2, y2-y1, lt = ',y1,y2, y2-y1, lt)
		#mapfull[y1:y2,:] = mapl;
		mapfull[y1:y1+lt,:] = mapl;
		y1 = y2;
	#--------------------------------
	o.map21 = mapfull;

#	print(mapfull)
#	print(np.shape(mapfull))
#	print(np.prod(np.shape(mapfull)))
	return 
#-----------------------------------------------
# create list of indices ranges that lead roughly to even load distribution in pool map.
def mkevenloadlist(trisize,Np):
	eloadl = [];
	dt=np.dtype('uint32');
	nblksize = np.zeros(Np,	dtype=dt);
	nntot = np.sum(trisize);
	dn = nntot/Np;
	#print(' dn  = ',dn)
	i = 0; i1 = 0; sm = 0; v2 = 1*dn; ichunk = 0;
	l = trisize.shape[0];
	while i < l-1 and ichunk < Np-1:
		sm += trisize[i];
		#print('i, v2 = ',i, v2)
		if sm >= v2:
			if i == i1:
				i2 = i+1;
				#print(' i=i1 ',i2)
			elif abs(sm-v2)<abs(sm-trisize[i]-v2):
				i2 = i;
				#print(' < ',i1)
			else:
				i2 = i-1;
				#print(' > ',i2)				
			eloadl.append([i1,i2]); ichunk += 1;
			#print('ichunk, dsize = ',ichunk,np.sum(trisize[i1:i2]))
			nblksize[ichunk-1] = np.sum(trisize[i1:i2]);
			i1 = i2; v2 += dn;	
			#i += 1;
		#else:
		i += 1;	
	eloadl.append([i1,l]); # last chunk
	#print('ichunk, dsize = ',ichunk+1,np.sum(trisize[i1:l]))
	#print(' size = ',np.sum(trisize))
	nblksize[Np-1] = np.sum(trisize[i1:l]);
	return eloadl, nblksize
#-----------------------------------------------


