
import globalvariables as o

from multiprocessing import Pool
from numpy import sqrt
from scipy.sparse import coo_matrix, diags


#****************************************************************
# sets global var: Hcsm, Hxsm, Hvsm, Hbsm, Hgsm, sft
#****************************************************************

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx, o.Np;
ld = o.ld
n1fsym,n1,n2,n3,ntot = o.n1fsym,o.n1,o.n2,o.n3,o.ntot;
listn2, listn3 = o.listn2, o.listn3
detuning = o.detuning;
corrtd = o.corrtd;

#****************************************************************
# flattens output of pool.map to make a list of elements
flat = lambda l: [item for sublist in l for item in sublist];
#****************************************************************
# calculates matrix elements of Hb for a given jj in range(0,n2)
# unexcited-sites-diagonal part:
def fHb(jj):
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
	for mj2 in range(m+1,mx+1):
		mj1 = mj2 - 1;	
		# transitions only for the same jj but adjuscent mi's.
		i1= n1fsym + (mj1-m-1)*n2 + jj;
		i2= n1fsym + (mj2-m-1)*n2 + jj;
		hbij = sqrt(mj2); 
		Hblp2ss.append([i1,i2,-ld*hbij]); 
	return Hblp2ss
# ------------------------------
# Photon: extra-basis: o.s
# special-sites-diagonal part: Pos = Photon other site
def fHbP2os(kchunk):
	hb2pos = [];
	for kk in range(kchunk[0],kchunk[1]):
		Pkk = o.Norm3l[kk];
		for mk2 in range(1,m+1):
			mk1 = mk2 - 1;
			jj1 = o.map32[kk,mk1];
			jj2 = o.map32[kk,mk2];
			Pjj1=o.Norm2l[jj1];
			Pjj2=o.Norm2l[jj2];
			for mj in range(m+1,mx+1):
				i1= n1fsym + (mj-m-1)*n2 + jj1; # coordinates
				i2= n1fsym + (mj-m-1)*n2 + jj2;
				yy = sqrt(mk2);
				xx = sqrt(Pjj1*Pjj2);
				xij=Pkk*(n-1)*yy/xx; # -ld displaced
				hb2pos.append([i1,i2,-ld*xij]);
	return hb2pos
#--------------------------------
# n=2 case: kchunk not available!
def fHbP2osn2(mj1): # pool map range(0,m), not m+1
	hb2posn2 = [];
	mj2 = mj1 + 1;
	for mj in range(m+1,mx+1):
		i1 = n1fsym + (mj-m-1)*n2 + mj1;
		i2 = n1fsym + (mj-m-1)*n2 + mj2;
		xij = sqrt(mj2);
		hb2posn2.append([i1,i2,-ld*xij]);
	return hb2posn2
# ------------------------------
# Photon: fsym_extra-basis (12):  s.s
def fHbP12ss(jj):
	# only mj1 = m, mj2 = m+1 because mj1 <= m & mj2 >= m+1
	mj1 = m; mj2 = m+1;	
	# transitions only for the same jj but adjuscent mi's.
	i1 = o.map21[jj,mj1];
	i2= n1fsym + jj;#i2= n1fsym + (mj2-m-1)*n2 + jj;
	Pii=o.Norm1l[i1];
	Pjj=o.Norm2l[jj];
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
				m11=o.Norm1l[i1];
				m12=o.Norm1l[i2];
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
				i1 = o.map32[kk,mj1];
				i2 = o.map32[kk,mj2];
				m11=o.Norm2l[i1];
				m12=o.Norm2l[i2];
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
			Nvj = o.Nv2l[jj] + mj
			i = n1fsym +(mj-m-1)*n2 + jj;
			Hvloc.append([i,i,Nvj])
	return Hvloc
#****************************************************************
# in exciton 0c block
# for j in range(0,n2):
def fHv0c(j):
	Hvloc = [];
	for mj in range(0,mx+1):
		Nvj = o.Nv2l[j] + mj
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
	print(' ====> calculating Hamiltonian ... ')
	pool=Pool(Np) # Np processes

	if ld>0:
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
				HbP2os=flat(HbP2os);
			else:
				HbP2os = [];
		elif (n>2):
			Hb2 = pool.map(fHb2,listn3);
			Hb2=flat(Hb2);
			if mx>m:
				HbP2os = pool.map(fHbP2os,listn3);
				HbP2os=flat(HbP2os);
			else:
				HbP2os = [];
		Hb1 = pool.map(fHb1,listn2);
		Hb1=flat(Hb1);
		if mx>m:
			HbP2ss = pool.map(fHbP2ss,range(0,n2));
			HbP2ss = flat(HbP2ss);
			HbP12ss = pool.map(fHbP12ss,range(0,n2));
		else:
			HbP12ss = []; HbP2ss=[];
		Hb2b = pool.map(fHb,range(0,n2));
		Hb2b=flat(Hb2b);
		Hb = Hb1+Hb2+Hb2b+HbP2os+HbP2ss+HbP12ss;
	else:
		# undisplaced basis: only electronically excited site is vib coupled:
		# calc of non-zero matrix elements of H_b = sum_n[c_n^deggar.c_n.b_n]
		Hb = pool.map(fHb,range(0,n2));
		Hb=flat(Hb);
		if 0:
			print("Hb = ", Hb )
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
	
	row=[];col=[];dat=[];
	for x in Hb:
		row.append(x[0])
		col.append(x[1])
		dat.append(x[2])
	Hbsm=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
	del Hb # free memory
	
	row=[];col=[];dat=[];
	for x in Hg:
		row.append(x[0])
		col.append(x[1])
		dat.append(x[2])
	Hgsm=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
	del Hg # free memory
	
	row=[];col=[];dat=[];
	for x in Hv:
		row.append(x[0])
		col.append(x[1])
		dat.append(x[2])
	o.Hvsm=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
	del Hv # free memory
	
	if (corrtd or detuning):
		row=[];col=[];dat=[];
		for x in Hc:
			row.append(x[0])
			col.append(x[1])
			dat.append(x[2])
		o.Hcsm=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
		del Hc # free memory
		row=[];col=[];dat=[];
		for x in Hx:
			row.append(x[0])
			col.append(x[1])
			dat.append(x[2])
		o.Hxsm=coo_matrix((dat, (row, col)), shape=(ntot, ntot))
		del Hx # free memory
		del row # free memory
		del col
		del dat
	# complete Hg & Hg matrices from upper triangular
	#Hv,Hc,Hx are diagonal so dont's do anything!
	o.Hbsm = Hbsm + Hbsm.transpose() - diags(Hbsm.diagonal(),0)
	o.Hgsm = Hgsm + Hgsm.transpose() - diags(Hgsm.diagonal(),0)


	# shifts due to half displaced basis
	# ld=1/2==> photon block: +1/4 *n*wv*lamb0**2 
	# ld=1/2==> exciton block: -3/4 *n*wv*lamb0**2 
	row=[];col=[];dat=[];
	disp1=n*ld**2;
	disp2=n*ld**2 - 2*ld;
	for i in range(n1):
		row.append(i)
		col.append(i)
		dat.append(disp1)
	for i in range(n2*(mx+1)):
		row.append(n1+i)
		col.append(n1+i)
		dat.append(disp2)
	o.sft=coo_matrix((dat, (row, col)), shape=(ntot, ntot));
	del row # free memory
	del col
	del dat
	# ------------------------------------------
	return














