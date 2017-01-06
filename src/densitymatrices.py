


import globalvariables as o
import numpy as np
from multiprocessing import Pool

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2 = o.n1,o.n2;
n1fsym = o.n1fsym;
nlmax = o.nlmax;
loopover= o.loopover;
lambin0 = o.lambin0;

listn2 = o.listn2;
listn3 = o.listn3;

# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

#****************************************************************
# symmetrise an upper or lower triangular dense matrix
def symmetrize(a):
    return a + a.T - np.diag(a.diagonal())
#****************************************************************

#****************************************************************
# Conditional density matrices
# dm0 ==> photon block;
# dm1 ==> exciton block,excited site
# dm2 ==> exciton block,unexcited site
#****************************************************************
# For dm0 & dm2 we need:
# mapmat0 with rows==>basjj & col==>mj; & matrix elem ==> basi
# norm_i for basi with elem ==> norm factors (# permutations)
# norm_jj for basjj & norm_k=norm_jjj for bask (N-2 sites)
#****************************************************************
# calc of dm0
#****************************************************************
def getrho0(jchunk):
	#dm0l = [[[0.0]*(m+1) for x in range(m+1)] for il in range(nlmax)];
	#dm0l = np.ndarray(dm0l);
	dm0l = np.zeros((nlmax,mx+1,mx+1));
	for jj in range(jchunk[0],jchunk[1]):
		m1 = o.Norm2l[jj];
		for mj1 in range(0,m+1):
			i1 = o.map21[jj,mj1];
			m11=o.Norm1l[i1];
			for mj2 in range(mj1,m+1): # upper triangular of rho0 only
				i2 = o.map21[jj,mj2];
				m12=o.Norm1l[i2];
				xx = np.sqrt(m11*m12);
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[i1][0]);
					cj=evec[i2][0];
					xij=m1*ci*cj/xx;
					dm0l[il,mj1,mj2] += xij;
	return dm0l
#****************************************************************


#--------------------------------------------
# btw extra and extra, type 1:
#--------------------------------------------
# diag in set jj: 
# rho[mj1,mj2]; mj1,mj2 special site occupation
def getrho0P2ss(jchunk):
	dm0l = np.zeros((nlmax,mx+1,mx+1));
	nm = n1fsym - (m+1)*n2;
	for jj in range(jchunk[0],jchunk[1]):
		nmj = nm + jj;
		for mj1 in range(m+1,mx+1):
			i1 = mj1*n2 + nmj;
			for mj2 in range(mj1,mx+1):
				i2 = mj1*n2 + nmj;
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[i1][0]);
					cj=evec[i2][0];
					dm0l[il,mj1,mj2] += ci*cj;
	return dm0l
#---------------------------
# btw extra and extra, type 2:
#----------------------------
# diag in special site mj: 
# rho[mk1,mk2]; mk1,mk2 in set jj
# that reduce to the same sym set kk (for n-2 sites).
def getrho0P2os(kchunk):
	dm0l = np.zeros((nlmax,mx+1,mx+1));
	nm = n1fsym - (m+1)*n2;
	for kk in range(kchunk[0],kchunk[1]):
		Pkkn = o.Norm3l[kk] * (n-1)/n;
		for mk1 in range(0,m+1):
			jj1 = o.map21[kk,mk1];
			Pj1 = o.Norm2l[jj1];
			nmj1 = nm + jj1;	
			for mk2 in range(mk1,m+1):
				jj2 = o.map21[kk,mk2];
				Pj2 = o.Norm2l[jj2];
				nmj2 = nm + jj2;	
				xij = Pkkn* np.sqrt(Pj1*Pj2);
				for mj in range(m+1,mx+1):
					i1 = mj*n2 + nmj1;
					i2 = mj*n2 + nmj2;
					for il, evalu, evec in o.eigvv:
						ci=np.conjugate(evec[i1][0]);
						cj=evec[i2][0];
						dm0l[il,mk1,mk2] += ci*cj*xij;
	return dm0l
#-------------------------
# n==2 version:
def getrho0P2osn2(mj1):
		dm0l = np.zeros((nlmax,mx+1,mx+1));
		nm = n1fsym - (m+1)*n2;
		xij0=0.5;#m1*(n-1)/(n*xx); m1=1,xx=1
		#for mj1 in range(0,m+1):
		nmj1 = nm + mj1;	
		for mj2 in range(mj1,m+1):
			nmj2 = nm + mj2;	
			for mj in range(m+1,mx+1):
				i1 = mj*n2 + nmj1;
				i2 = mj*n2 + nmj2;
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[i1][0]);
					cj=evec[i2][0];
					dm0l[il,mk1,mk2] += ci*cj*xij0;
		return dm0l
#--------------------------------------------
# btw fsym and extra:
#--------------------------------------------
# special site mj > m so contrib to only rho[mj1,mj],
# with mj1 in fsym set i that reduce to the same sym set
# jj as for extra state's sym set part.
def getrho0P2fe(jchunk):
	dm0l = np.zeros((nlmax,mx+1,mx+1));
	nm = n1fsym - (m+1)*n2;
	for jj in range(jchunk[0],jchunk[1]):
		Pjj = o.Norm2l[jj];
		nmj = nm + jj;
		for mj1 in range(0,m+1):
			i1 = o.map21[jj,mj1];
			Pi1 = o.Norm1l[i1];
			for mj2 in range(m+1,mx+1): # upper triangular of rho0 only
				i2 = nmj + mj2*n2;
				xij = np.sqrt(n*Pjj/Pi1);
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[i1][0]);
					cj=evec[i2][0];
					dm0l[il,mj1,mj2] += ci*cj*xij;
	return dm0l

#****************************************************************
# calc of dm1
#****************************************************************
def getrho1(jchunk):
	dm1l = np.zeros((nlmax,mx+1,mx+1));
	for i in range(jchunk[0],jchunk[1]):
		for mi in range(0,mx+1):
			ii = mi*n2 + i # index of basis state in this block
			for mj in range(mi,mx+1): # only upper triangular
				jj = mj*n2 + i  # only for j = i terms contrib
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[n1+ii][0]);
					cj=evec[n1+jj][0]
					xij = ci*cj/n;
					dm1l[il,mi,mj] += xij
	return dm1l
#****************************************************************
# calc of dm2
#****************************************************************
def getrho2(kchunk):
	dm2l = np.zeros((nlmax,m+1,m+1));
	for kk in range(kchunk[0],kchunk[1]):
		m1n = o.Norm3l[kk] * (n-1)/n;
		for mj1 in range(0,m+1):
			i1 = o.map21[kk,mj1];
			m11=o.Norm2l[i1];
			for mj2 in range(mj1,m+1): # upper triangular of rho0 only
				i2 = o.map21[kk,mj2];
				m12=o.Norm2l[i2];
				xx = np.sqrt(m11*m12);
				xij0=m1n/xx;
				for il, evalu, evec in o.eigvv:
					for m2 in range(mx+1):
						ii1=m2*n2 + i1; # diag in m2 contrib only
						ii2=m2*n2 + i2;
						ci=np.conjugate(evec[n1+ii1][0]);
						cj=evec[n1+ii2][0];
						xij = xij0*ci*cj;
						dm2l[il,mj1,mj2] += xij
	return dm2l
#****************************************************************
# dm2 for n=2 case
def getrho2n2():
	dm2l = np.zeros((nlmax,m+1,m+1));
	xij0=0.5;#m1*(n-1)/(n*xx); m1=1,xx=1
	for mj1 in range(0,m+1):
			for mj2 in range(mj1,m+1): # upper triangular of rho0 only
				for il, evalu, evec in o.eigvv:
					for m2 in range(mx+1):
						ii1=m2*n2 + mj1; # diag in m2 contrib only
						ii2=m2*n2 + mj2; 
						ci=np.conjugate(evec[n1+ii1][0]);
						cj=evec[n1+ii2][0];
						xij = xij0*ci*cj;
						dm2l[il,mj1,mj2] += xij
	return dm2l
#****************************************************************















def cdms():
	print(' ====> calculating conditional vib dms... ')
	# -----------------------------
	if loopover == 'lambda0':
	# parameters:
		# param = np.array([n,m,wv,wr,lamb]);
		param = np.array([n,m,wv,wr])
	# descriptive output file names
		fout0=dumy+"/rho0_vib-N="+str(n)+"-M="+str(m)+"-wr="+str('%.1f'%(wr))+".dat"
		fout1=dumy+"/rho1_vib-N="+str(n)+"-M="+str(m)+"-wr="+str('%.1f'%(wr))+".dat"
		fout2=dumy+"/rho2_vib-N="+str(n)+"-M="+str(m)+"-wr="+str('%.1f'%(wr))+".dat"
	else: # if loopover == 'wr':
	# parameters:
		# param = np.array([n,m,wv,wr,lamb]);
		param = np.array([n,m,wv,lambda0])
	# descriptive output file names
		fout0=dumy+"/rho0_vib-N="+str(n)+"-M="+str(m)+"-l0="+str('%.1f'%(lambda0))+".dat"
		fout1=dumy+"/rho1_vib-N="+str(n)+"-M="+str(m)+"-l0="+str('%.1f'%(lambda0))+".dat"
		fout2=dumy+"/rho2_vib-N="+str(n)+"-M="+str(m)+"-l0="+str('%.1f'%(lambda0))+".dat"
	# -----------------------------

	#***************************************************
	print('      calculating conditional reduced density matrices ')
	# -----------------------------
	# n=1 case:
	# -----------------------------
	if n==1:
		dm0, dm1 = cdmsn1();
		writedms(dm0,param,lambin0,fout0);
		writedms(dm1,param,lambin0,fout1);
		return
	# -----------------------------	
	# n > 1 cass:
	# -----------------------------
	pool=Pool(Np);	
	#*********************************
	# calculate the reduced density matrix
	#  of vib DOF of a single site
	#******************************
	# calc of dm0
	#*******************************
	print('       calculating dm0 ... ')
	# original fully sym basis, now first n1fsym basis in photon block
	dm0s = pool.map(getrho0,listn2)  # note list2
	dm0 = sum(dm0s,0) # sum along 0 ==> pool map
	# due to extra basis ---------------
	# extra-extra: spec site (type 1)
	dm0s = pool.map(getrho0P2ss,listn2);
	dm0 += sum(dm0s,0);
	# extra-extra: other site (type 2)
	if n==2:
		dm0s = pool.map(getrho0P2osn2,range(0,m+1))
		dm0 += sum(dm0s,0);
	elif n>2:
		dm0s = pool.map(getrho0P2os,listn3)
		dm0 += sum(dm0s,0);
	# fsym-extra:
	dm0s = pool.map(getrho0P2fe,listn2)
	dm0 += sum(dm0s,0);
	del dm0s # free memory
	# -----------------------------------
	writedms(dm0,param,lambin0,fout0);
	#*********************************
	# calc of dm0
	#*********************************
	# calc dm1: layers <==> processes
	print('       calculating dm1 ... ')
	dm1s = pool.map(getrho1,listn2)  # note list2 again
	dm1 = sum(dm1s,0) # sum along which dim, 3? or 0?
	del dm1s # free memory
	writedms(dm1,param,lambin0,fout1);
	#***********************************
	# calc dm2: layers <==> processes
	# uses listn2
	print('       calculating dm2 ... ')
	if n>2:
		dm2s = pool.map(getrho2,listn3)  # note list3!
		dm2 = sum(dm2s,0) # sum along dim=0
		del dm2s # free memory
	else: # n = 2 case
		dm2 = getrho2n2()
	writedms(dm2,param,lambin0,fout2);	
	pool	.close();
	return
#-------------------------------------
# for n=1 
#-------------------------------------
def cdmsn1():
	#-------------------------------------
	def getrhon1(p):
		m = mx; # for n=1, vib cutoff = mx
		dm1l = np.zeros((nlmax,m+1,m+1));
		for mi in range(0,m+1):
			for mj in range(mi,m+1): # only upper triangular
				for il, evalu, evec in o.eigvv:
					ci=np.conjugate(evec[p+mi][0]);
					cj=evec[p+mj][0]
					xij = ci*cj;
					dm1l[il,mi,mj] += xij
		return dm1l
	#-------------------------------------
	dm0 = getrhon1(0);
	dm1 = getrhon1(n1);
	return dm0, dm1
#-------------------------------------
# for writing output files 
#-------------------------------------
def writedms(dm,param,lambin0,fout):
	f=open(fout,'ab')
	np.savetxt(f,param[None],  fmt='%5i, %5i, %10.6f, %10.6f')
	np.savetxt(f,lambin0[None],  fmt='%10.6f',delimiter=',')
	for x in dm:
		y=symmetrize(x)
		np.savetxt(f,y, fmt='%15.10f', delimiter=',')
	f.close()
	return
#-------------------------------------

	
