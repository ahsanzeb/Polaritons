
import globalvariables as o

from multiprocessing import Pool
from numpy import sqrt
from scipy.sparse import coo_matrix, diags
import numpy as np
import gc
from memtime import memtime

#****************************************************************
# sets global var: Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
#****************************************************************

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx, o.Np;
ld = o.ld;
# ld  = 0.5;
n1fsym,n1,n2,n3,ntot = o.n1fsym,o.n1,o.n2,o.n3,o.ntot;
listn1fsym, listn2, listn3 = o.listn1fsym, o.listn2, o.listn3
detuning = o.detuning;
corrtd = o.corrtd;

#print('hamiltonian: ld = ',ld)

# print('mod ham: n1,n2,ntot, mx = ',n1,n2,ntot, mx )
#print("n,m,mx =",n,m,mx )
#****************************************************************
# flattens output of pool.map to make a list of elements
flat = lambda l: [item for sublist in l for item in sublist];
#****************************************************************
# calculates matrix elements of Hb for a given jj in range(0,n2)
# unexcited-sites-diagonal part:
def fHb(jj):
# 	print('fHb: n1,n2,ntot, mx = ',n1,n2,ntot, mx)
	Hbloc=np.zeros((3,mx+1));
	for mi in range(1,mx+1):
		# transitions only for the same jj but adjuscent mi's.
		i=(mi-1)*n2+jj # coor in lower right block
		j=mi*n2 + jj 
		hbij = sqrt(mi); 
		#Hbloc.append([n1+i,n1+j,(1-ld)*hbij]);
		Hbloc[:,mi]= n1+i,n1+j,(1-ld)*hbij
	return Hbloc

# ------------------------------
# 16-12-16 for absorption extra basis
# 'not-special'-sites-diagonal part: Pss = Photon special site
# Photon: extra-basis: s.s
def fHbP2ss(jj):
	Hblp2ss=np.zeros((3, mx-m))
	mn2j =n1fsym -(m+1)*n2 + jj;
	ind = 0;
	for mj2 in range(m+1,mx+1):
		mj1 = mj2 - 1;	
		# transitions only for the same jj but adjuscent mi's.
		i1= mj1*n2 +mn2j;
		i2= mj2*n2 +mn2j;
		hbij = sqrt(mj2); 
		Hblp2ss[:,ind] = i1,i2,-ld*hbij
		ind += 1;
	return Hblp2ss
# ------------------------------
# Photon: extra-basis: o.s
# special-sites-diagonal part: Pos = Photon other site
def fHbP2os(kchunk):
	hb2pos = np.zeros((3,m*(mx-m)*(kchunk[1]-kchunk[0]) ))
	mn2 = n1fsym -(m+1)*n2;
	ind = 0;
	for kk in range(kchunk[0],kchunk[1]):
		Pkk = o.Norm3l[kk];
		for mk2 in range(1,m+1):
			mk1 = mk2 - 1;
			jj1 = o.map21[kk,mk1];
			jj2 = o.map21[kk,mk2];
			Pjj1= o.Norm2l[jj1];
			Pjj2= o.Norm2l[jj2];
			yy = sqrt(mk2);
			xx = sqrt(Pjj1*Pjj2);
			xij=Pkk*(n-1)*yy/xx; # -ld displaced
			mnj1 = mn2 + jj1;
			mnj2 = mn2 + jj2;
			for mj in range(m+1,mx+1):
				i1= mj*n2 + mnj1;
				i2= mj*n2 + mnj2;
				hb2pos[:,ind] = i1,i2,-ld*xij
				ind += 1;
	return hb2pos
#--------------------------------
# n=2 case: kchunk not available!
def fHbP2osn2(mj1): # pool map range(0,m), not m+1
	hb2posn2 = np.zeros((3,mx-m));
	ind = 0;
	mj2 = mj1 + 1;
	mn2j1 = n1fsym-(m+1)*n2 + mj1
	mn2j2 = n1fsym-(m+1)*n2 + mj2
	xij = sqrt(mj2);
	for mj in range(m+1,mx+1):
		i1 = mj*n2 + mn2j1;
		i2 = mj*n2 + mn2j2;
		hb2posn2[:,ind] = i1,i2,-ld*xij
		ind += 1;
	return hb2posn2
# ------------------------------
# Photon: fsym_extra-basis (12):  s.s
def fHbP12ss(jj):
	x = np.zeros((3,1));
	# only mj1 = m, mj2 = m+1 because mj1 <= m & mj2 >= m+1
	mj1 = m; mj2 = m+1;	
	# transitions only for the same jj but adjuscent mi's.
	i1 = o.map21[jj,mj1];
	i2 = n1fsym + jj;#i2= n1fsym + (mj2-m-1)*n2 + jj;
	Pii= o.Norm1l[i1];
	Pjj= o.Norm2l[jj];
	hbij = sqrt(mj2*n*Pjj/Pii);
	x[:,0] = i1,i2,-ld*hbij 
	return x
# ------------------------------
# ------------------------------

# In half displaced basis hdb (displaced by ld),
# H_{\lambda} has matrix elements in both
# 		photon and exciton blocks
# 	H ===> H 	+	N*wv*l0^2 
#						- 2wv*l0*ld*Sum_i[c_i^dag*ci] 
#						-	wv*ld*Sum_i[bi^dag + bi]
# ------------------------------
# photon block:
def fHb1(jchunk):
	hb1=np.zeros((3,m*(jchunk[1]-jchunk[0]) ));
	ind = 0;
	for jj in range(jchunk[0],jchunk[1]):
		m1 = o.Norm2l[jj];
		for mj1 in range(0,m):
				mj2 = mj1 + 1;
				i1 = o.map21[jj,mj1];
				i2 = o.map21[jj,mj2];
				m11= o.Norm1l[i1];
				m12= o.Norm1l[i2];
				yy = sqrt(mj2);
				xx = sqrt(m11*m12);
				xij = -ld*m1*n*yy/xx; # -ld displaced 
				hb1[:,ind]= i1,i2,xij
				ind += 1;
	return hb1
# ------------------------------
# exciton block: excited-site-diagonal part
# ------------------------------
def fHb2(kchunk):
	hb2=np.zeros((3,m*(mx+1)*(kchunk[1]-kchunk[0]) ));
	ind = 0;
	for kk in range(kchunk[0],kchunk[1]):
		m1 = o.Norm3l[kk];
		for mj1 in range(0,m):
				mj2 = mj1 +1;
				i1 = o.map21[kk,mj1];
				i2 = o.map21[kk,mj2];
				m11= o.Norm2l[i1];
				m12= o.Norm2l[i2];
				yy = sqrt(mj2);
				xx = sqrt(m11*m12);
				xij=-ld*m1*(n-1)*yy/xx; # -ld displaced
				for m2 in range(0,mx+1):
					ii1=m2*n2 + i1; # diag in m2 contrib only
					ii2=m2*n2 + i2;
					hb2[:,ind]= n1+ii1,n1+ii2,xij
					ind += 1;
	return hb2
# ------------------------------	
def fHb2n2(mj1):
	hb2 = np.zeros((3,mx+1))
	mj2 = mj1+1;
	i1 = mj1;
	i2 = mj2;
	yy = sqrt(mj2);
	xij=-ld*yy; # -ld displaced
	for m2 in range(0,mx+1):
		ii1=m2*n2 + i1; # diag in m2 contrib only
		ii2=m2*n2 + i2;
		hb2[:,m2] = n1+ii1,n1+ii2,xij
	return hb2
# ------------------------------	

#****************************************************************
# calculate matrix elements of Hg using mapping21 (nparray)
def getHg(jchunk):
	Hgloc = np.zeros((3, (mx+1)*(jchunk[1]-jchunk[0])));
	ind = 0;
	for jj in range(jchunk[0],jchunk[1]):
		m2 = o.Norm2l[jj];
		for mj in range(0,m+1):
			j= n1 + mj*n2 + jj # coor in uper right block
			ii = o.map21[jj,mj]
			m1 = o.Norm1l[ii];
			hgij=sqrt(n*m2/m1)
			Hgloc[:,ind] = ii,j,hgij
			ind += 1;
		# extra basis: diag in mj,jj
		for mj in range(m+1,mx+1): 
			# coor in exciton block
			j = n1 + mj*n2 + jj; 
			# shift j by -(n1+n2*(m+1)) to get relative coor
			# then shift by + n1fsym
			ii = n1fsym +(mj-m-1)*n2 + jj;
			Hgloc[:,ind] = ii,j,1
			ind += 1;
	return Hgloc
#****************************************************************
# in photon 1c block
# for i in range(0,n1fsym):
def fHv1c(ichunk):
	Hvloc = np.zeros((3,ichunk[1]-ichunk[0]))
	ind=0;
	for i in range(ichunk[0],ichunk[1]):
		Nvi = o.Nv1l[i];
		Hvloc[:,ind] = i,i,Nvi
		ind += 1;
	return Hvloc

def fHv1cnew():
	Hv1c = np.zeros((3,n1fsym));
	Hv1c[0,:] = range(n1fsym)
	Hv1c[1,:] = range(n1fsym)
	Hv1c[2,:] = o.Nv1l[:];
	return Hv1c

# Hv for extra states in photon block:
def fHv1cExtra(jchunk):
	Hvloc = np.zeros((3,(mx-m)*(jchunk[1]-jchunk[0])))
	ind=0;
	for jj in range(jchunk[0],jchunk[1]):
		for mj in range(m+1,mx+1):
			Nvj = o.Nv1l[jj] + mj
			i = n1fsym +(mj-m-1)*n2 + jj;
			Hvloc[:,ind] = i,i,Nvj
			ind += 1;
	return Hvloc
#****************************************************************
# in exciton 0c block
def fHv0c(jchunk):
	Hvloc = np.zeros((3,(mx+1)*(jchunk[1]-jchunk[0])))
	ind=0;
	for j in range(jchunk[0],jchunk[1]):
		for mj in range(0,mx+1):
			Nvj = o.Nv1l[j] + mj
			jj = n1 + mj*n2 + j
			Hvloc[:,ind] = jj,jj,Nvj
			ind += 1;
	return Hvloc
#****************************************************************
# Hc & Hx
#for i in range(0,n1):
def fHc(i):
	return [i,i,1]
#****************************************************************
#Hx = [];
#for i in range(n1,ntot):
def fHx(i):
	return [i,i,1]
#****************************************************************

def fHcnew():
	Hc = np.zeros((3,n1));
	Hc[0,:] = range(n1)
	Hc[1,:] = range(n1)
	Hc[2,:] = 1;
	return Hc

def fHxnew():
	Hx = np.zeros((3,ntot-n1));
	Hx[0,:] = range(n1,ntot)
	Hx[1,:] = range(n1,ntot)
	Hx[2,:] = 1;
	return Hx

def fsftm():
	# shifts due to half displaced basis
	disp1=n*ld**2; # photon block
	disp2=n*ld**2 - 2*ld; # exciton block
	# not -2*n*ld: sum_i[cidag*ci (-2lam)] = -2lam*sum_i[cidag*ci] = -2lam
	sftm = np.zeros((3,ntot));
	sftm[0,:] = range(ntot);
	sftm[1,:] = range(ntot);
	sftm[2,:n1] = disp1;	
	sftm[2,n1:] = disp2;	
	return sftm
		



#****************************************************************
#  Hamiltonain: parallel
#****************************************************************
def hamilt():

	#print('listn2 = ',listn2)
	#print('listn3 = ',listn3)	#****************************************************************
	# calculate Hamiltonain and set global variables:
	# Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
	#****************************************************************
	print(' ====> calculating Hamiltonian ... ');
# -------------------------
#  n =1 case:
# -------------------------
	if n == 1:
		hamiltn1(o.mx);
		return
# -------------------------
# n > 1 cases:
# -------------------------
	pool=Pool(Np) # Np processes
	# -----------------------------------------
	if ld>0:
		# -----------------------------------------
		# displaced basis:
		# H_{\lambda} has matrix elements in both photon and exciton blocks
		# 	H ===> H 	+	N*wv*l0^2 
		#						- 2wv*l0*ld*Sum_i[c_i^dag*ci] 
		#						-	wv*ld*Sum_i[bi^dag + bi]
		# computing matrix elements: 
		# 	extra complication due to photon block's extra basis states
		if(n==2):
			Hb2 = pool.map(fHb2n2,range(m)); # not range(m+1), it's row indeces for b
			if mx>m:
				HbP2os = pool.map(fHbP2osn2,range(m));
			else:
				HbP2os = [];
		elif (n>2):
			Hb2 = pool.map(fHb2,listn3);
			if mx>m:
				HbP2os = pool.map(fHbP2os,listn3);
			else:
				HbP2os = [];

		Hb1 = pool.map(fHb1,listn2);
		if mx>m:
			HbP2ss = pool.map(fHbP2ss,range(0,n2));
			HbP12ss = pool.map(fHbP12ss,range(0,n2));
		else:
			HbP12ss = []; HbP2ss=[];
		Hb2b = pool.map(fHb,range(0,n2));
		Hb = Hb1+Hb2+Hb2b+HbP2os+HbP2ss+HbP12ss;
	# -----------------------------------------
	else:
		Hb = [];
		# undisplaced basis: only electronically excited site is vib coupled:
		# calc of non-zero matrix elements of H_b = sum_n[c_n^deggar.c_n.b_n]
		Hb = pool.map(fHb,range(0,n2));

	Hb = np.hstack(Hb);
	o.Hbsm = coomatnp(Hb,1); 
	del Hb; gc.collect()
	memtime('Hb');
	print("        Hb calculated...")

	Hg=pool.map(getHg,listn2);
	Hg = np.hstack(Hg);
	o.Hgsm = coomatnp(Hg,1); del Hg; gc.collect()
	memtime('Hg');
	print("        Hg calculated...")

	# in photon 1c block
	#Hv1c= pool.map(fHv1c,listn1fsym); # range(0,n1fsym)
	Hv1c= [fHv1cnew()];
	Hv1cExtra= pool.map(fHv1cExtra,listn2)
	# in exciton 0c block
	Hv0c=pool.map(fHv0c,listn2);
	Hv = Hv1c + Hv1cExtra + Hv0c;
	Hv = np.hstack(Hv);
	o.Hvsm = coomatnp(Hv,0); del Hv; gc.collect()
	memtime('Hv');
	print("        Hv calculated...")
	
	# Hc & Hx
	if (corrtd or detuning): # corrtd or 
		Hc= fHcnew()
		o.Hcsm = coomatnp(Hc,0); del Hc; gc.collect()
		Hx= fHxnew();		
		o.Hxsm = coomatnp(Hx,0); del Hx; gc.collect()
	else:
		o.Hcsm= coo_matrix(([], ([], [])), shape=(ntot, ntot))	
		o.Hxsm= coo_matrix(([], ([], [])), shape=(ntot, ntot))	
	memtime('Hcx');
	print("        Hc,x calculated... ")

	pool.close() # this pool does not know various variables calculated after its start

	if abs(ld) > 1e-5:
		sftm = fsftm();
		o.sft= coomatnp(sftm,0); del sftm; gc.collect()
	else:
		o.sft= coo_matrix(([], ([], [])), shape=(ntot, ntot))
	memtime('sft');
	# ------------------------------------------
	return



def hamiltn1(mx):
	# ------------------------
	# undisplaced phonon basis will be used for n=1
	# ------------------------
	Hg = [];
	for i in range(mx+1):
		Hg.append([i,n1+i,1]);
	o.Hgsm = coomat(Hg,1); del Hg;
	Hb = [];
	for i in range(1,mx+1):
		Hb.append([i-1,i,sqrt(i)]);
	o.Hbsm = coomat(Hb,1);  del Hb;
	Hv = [];
	for i in range(0,mx+1):
		Hv.append([i,i,i]);
		Hv.append([n1+i,n1+i,i]);
	o.Hvsm = coomat(Hv,0);  del Hv;
	if (corrtd or detuning): #corrtd or 
		Hc = [];
		for i in range(0,mx+1):
			Hc.append([i,i,1]);
		o.Hcsm = coomat(Hc,0);  del Hc;
		Hx = [];
		for i in range(0,mx+1):
			Hx.append([n1+i,n1+i,1]);
		o.Hxsm = coomat(Hx,0);  del Hx;
	o.sft= coomat([],0); # undisplaced basis
	# ------------------------
	return

def coomatnp(x,sym):
	Y = coo_matrix((x[2,:], (x[0,:], x[1,:])), shape=(ntot, ntot))
	if sym:
		# complete matrices from upper/lower triangular
		Y = Y + Y.transpose() - diags(Y.diagonal(),0)
	return Y

def coomat(X,sym):
		row=[];col=[];dat=[];
		for x in X:
			row.append(x[0])
			col.append(x[1])
			dat.append(x[2])
		Y=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
		if sym:
			# complete matrices from upper/lower triangular
			Y = Y + Y.transpose() - diags(Y.diagonal(),0)
		return Y

