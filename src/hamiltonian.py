
import globalvariables as o

from multiprocessing import Pool
from numpy import sqrt
from scipy.sparse import coo_matrix, diags
import numpy as np

#****************************************************************
# sets global var: Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
#****************************************************************

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx, o.Np;
ld = o.ld;
# ld  = 0.5;
n1fsym,n1,n2,n3,ntot = o.n1fsym,o.n1,o.n2,o.n3,o.ntot;
listn2, listn3 = o.listn2, o.listn3
detuning = o.detuning;
corrtd = o.corrtd;

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
	Hbloc=[];
	for mi in range(1,mx+1):
		# transitions only for the same jj but adjuscent mi's.
		i=(mi-1)*n2+jj # coor in lower right block
		j=mi*n2 + jj 
		hbij = sqrt(mi); 
		Hbloc.append([n1+i,n1+j,(1-ld)*hbij]);
	return Hbloc

# ------------------------------
# 16-12-16 for absorption extra basis
# 'not-special'-sites-diagonal part: Pss = Photon special site
# Photon: extra-basis: s.s
def fHbP2ss(jj):
	Hblp2ss=[];
	mn2j =n1fsym -(m+1)*n2 + jj;
	for mj2 in range(m+1,mx+1):
		mj1 = mj2 - 1;	
		# transitions only for the same jj but adjuscent mi's.
		i1= mj1*n2 +mn2j;
		i2= mj2*n2 +mn2j;
		hbij = sqrt(mj2); 
		Hblp2ss.append([i1,i2,-ld*hbij]); 
	return Hblp2ss
# ------------------------------
# Photon: extra-basis: o.s
# special-sites-diagonal part: Pos = Photon other site
def fHbP2os(kchunk):
	hb2pos = [];
	mn2 = n1fsym -(m+1)*n2
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
				hb2pos.append([i1,i2,-ld*xij]);
	return hb2pos
#--------------------------------
# n=2 case: kchunk not available!
def fHbP2osn2(mj1): # pool map range(0,m), not m+1
	hb2posn2 = [];
	mj2 = mj1 + 1;
	mn2j1 = n1fsym-(m+1)*n2 + mj1
	mn2j2 = n1fsym-(m+1)*n2 + mj2
	xij = sqrt(mj2);
	for mj in range(m+1,mx+1):
		i1 = mj*n2 + mn2j1;
		i2 = mj*n2 + mn2j2;
		hb2posn2.append([i1,i2,-ld*xij]);
	return hb2posn2
# ------------------------------
# Photon: fsym_extra-basis (12):  s.s
def fHbP12ss(jj):
	# only mj1 = m, mj2 = m+1 because mj1 <= m & mj2 >= m+1
	mj1 = m; mj2 = m+1;	
	# transitions only for the same jj but adjuscent mi's.
	i1 = o.map21[jj,mj1];
	i2 = n1fsym + jj;#i2= n1fsym + (mj2-m-1)*n2 + jj;
	Pii= o.Norm1l[i1];
	Pjj= o.Norm2l[jj];
	hbij = sqrt(mj2*n*Pjj/Pii); 
	return [i1,i2,-ld*hbij];
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
	hb1 = [];
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
				hb1.append([i1,i2,xij]);
	return hb1
# ------------------------------
# exciton block: excited-site-diagonal part
# ------------------------------
def fHb2(kchunk):
	hb2 = [];
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
					hb2.append([n1+ii1,n1+ii2,xij]);
	return hb2
# ------------------------------	
def fHb2n2(mj1):
	hb2 = [];
	mj2 = mj1+1;
	i1 = mj1;
	i2 = mj2;
	yy = sqrt(mj2);
	xij=-ld*yy; # -ld displaced
	for m2 in range(0,mx+1):
		ii1=m2*n2 + i1; # diag in m2 contrib only
		ii2=m2*n2 + i2;
		hb2.append([n1+ii1,n1+ii2,xij]);
	return hb2
# ------------------------------	

#****************************************************************
# calculate matrix elements of Hg using mapping21 (nparray)
def getHg(jchunk):
	Hgloc = [];
	for jj in range(jchunk[0],jchunk[1]):
		for mj in range(0,m+1):
			j= n1 + mj*n2 + jj # coor in uper right block
			ii = o.map21[jj,mj]
			m1 = o.Norm1l[ii];
			m2 = o.Norm2l[jj];
			hgij=sqrt(n*m2/m1)
			Hgloc.append([ii,j,hgij]);
		# extra basis: diag in mj,jj
		for mj in range(m+1,mx+1): 
			# coor in exciton block
			j = n1 + mj*n2 + jj; 
			# shift j by -(n1+n2*(m+1)) to get relative coor
			# then shift by + n1fsym
			ii = n1fsym +(mj-m-1)*n2 + jj;
			Hgloc.append([ii,j,1]); 
	return Hgloc
#****************************************************************
# in photon 1c block
# for i in range(0,n1fsym):
def fHv1c(i): 
	Nvi = o.Nv1l[i]
	return [i,i,Nvi]	
# Hv for extra states in photon block:
def fHv1cExtra(jchunk):
	Hvloc = [];
	for jj in range(jchunk[0],jchunk[1]):
		for mj in range(m+1,mx+1):
			Nvj = o.Nv1l[jj] + mj
			i = n1fsym +(mj-m-1)*n2 + jj;
			Hvloc.append([i,i,Nvj])
	return Hvloc
#****************************************************************
# in exciton 0c block
# for j in range(0,n2):
def fHv0c(j):
	Hvloc = [];
	for mj in range(0,mx+1):
		Nvj = o.Nv1l[j] + mj
		jj = n1 + mj*n2 + j
		Hvloc.append([jj,jj,Nvj])
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




#****************************************************************
#  Hamiltonain: parallel
#****************************************************************
def hamilt():
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
			Hb2=flat(Hb2);
			if mx>m:
				HbP2os = pool.map(fHbP2osn2,range(m));
				HbP2os = flat(HbP2os);
			else:
				HbP2os = [];
		elif (n>2):
			Hb2 = pool.map(fHb2,listn3);
			Hb2=flat(Hb2);
			if mx>m:
				HbP2os = pool.map(fHbP2os,listn3);
				HbP2os = flat(HbP2os);
			else:
				HbP2os = [];

		Hb1 = pool.map(fHb1,listn2);
		Hb1 = flat(Hb1);
		if mx>m:
			HbP2ss = pool.map(fHbP2ss,range(0,n2));
			HbP2ss = flat(HbP2ss);
			HbP12ss = pool.map(fHbP12ss,range(0,n2));
		else:
			HbP12ss = []; HbP2ss=[];

		Hb2b = pool.map(fHb,range(0,n2));
		Hb2b=flat(Hb2b);
		Hb = Hb1+Hb2+Hb2b+HbP2os+HbP2ss+HbP12ss;
		if 0:
			print('------------');
			print(Hb1)
			print(Hb2)
			print(Hb2b)
			print(HbP2os)
			print(HbP2ss)
			print(HbP12ss)
	# -----------------------------------------
	else:
		Hb = [];
		# undisplaced basis: only electronically excited site is vib coupled:
		# calc of non-zero matrix elements of H_b = sum_n[c_n^deggar.c_n.b_n]
		Hb = pool.map(fHb,range(0,n2));
		Hb=flat(Hb);
	# -----------------------------------------
	print("        Hb calculated...")

	Hg=pool.map(getHg,listn2);
	Hg=flat(Hg)
	print("        Hg calculated...")

	# in photon 1c block
	Hv1c= pool.map(fHv1c,range(0,n1fsym))
	Hv1cExtra= pool.map(fHv1cExtra,listn2)
	# in exciton 0c block
	Hv0c=pool.map(fHv0c,range(0,n2));
	# flat and join lists
	Hv = Hv1c + flat(Hv1cExtra) + flat(Hv0c)
	print("        Hv calculated...")
	
	# Hc & Hx
	if (corrtd or detuning):
		Hc= pool.map(fHc,range(0,n1))
		Hx= pool.map(fHx,range(n1,ntot))
		print("        Hc,x calculated... ")
	
	pool.close() # this pool does not know various variables calculated after its start
	#****************************************************************
	# form complete Hamiltonain and set global variables:
	# Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
	#****************************************************************
	o.Hgsm = coomat(Hg,1); del Hg;
	o.Hbsm = coomat(Hb,1); del Hb;
	o.Hvsm = coomat(Hv,0); del Hv;
	if (corrtd or detuning):
		o.Hcsm = coomat(Hc,0); del Hc
		o.Hxsm = coomat(Hx,0); del Hx
	# shifts due to half displaced basis
	disp1=n*ld**2; # photon block
	disp2=n*ld**2 - 2*ld; # exciton block
	# not -2*n*ld: sum_i[cidag*ci (-2lam)] = -2lam*sum_i[cidag*ci] = -2lam
	sftm = [];
	for i in range(n1):
		sftm.append([i,i,disp1])
	for i in range(n1,ntot):
		sftm.append([i,i,disp2])
	o.sft= coomat(sftm,0); del sftm;
	# ------------------------------------------
	#print('o.sft='); print(o.sft)
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
	if (corrtd or detuning):
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


