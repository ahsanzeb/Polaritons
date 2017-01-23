# kappa, gamma multiplication
#when getting Gamma, use appropriate weights, Nv1l

import globalvariables as o
import numpy as np;
from math import factorial
from scipy.integrate import ode
from numpy import pi, exp
from scipy.sparse import coo_matrix
from scipy import linalg
import matplotlib.pyplot as plt
from multiprocessing import Pool
import decimal

from memtime import memtime

# define local to avoid writing o. every time!
# at risk of setting these at the time of first import of this module only.
##### any update in the values of these variables will not be considered by the functions of this module, because they use the first-import-time's set values.
########### consider this if loops over lambda/wr/n for absorption calculations is tried.
usenumpy = o.usenumpy;
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2,ntot = o.n1,o.n2,o.ntot;
n1fsym = o.n1fsym;
ld = o.ld;
tolr, itermax = o.tolr, o.itermax 
printstep = o.printstep;

gamma, kappa =	 o.gamma, o.kappa
if abs(gamma-kappa) > 1e-5: dkapa = 1;
else: dkapa = 0;

e1,e2,dt,tf = 	o.e1, o.e2, o.dt, o.tf
nwmax,ntmax =	o.nwmax, o.ntmax
show = 	o.show

kl= kappa/2;

nlmax = o.nlmax;
loopover= o.loopover;
lambda0 = o.lambda0;
lambin0 =	o.lambin0; # list of lam or wr values


# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

#----------------------------------
#	absorption using td: functions
	#	steps:
		# start with |P(t=0)> = a^dag|0x,0p,0v>, 
		# evolve to get |P(t)> & get the correlation function <P(0)|P(t)>
		# Fourier transform to get the absorption
#--------------------------
def gettd(il):
	# time evolves the state with only a photon present
	# psi0 = |P(t=0)> = a^dag|0x,0p,0v> = [1,0,0,......,0];
	#--------------------------		
	# making internal function
	# to avoid making list of H for loopover lam/wr
	#--------------------------
	kappa2 = kappa/2.0; gamma2 = gamma/2.0;
	def fInteg(t, y):
		# gives RHS of td-schrodinger equation with losses:
		# gives [-iota*H - (kappa*PhotProjector + gamma*ExcProjector)]psi(t)
		if dkapa: hdecay = np.concatenate( (kappa2*y[range(n1)],gamma2*y[range(n1,ntot)]) );
		else: hdecay = kappa2*y;	
		Hpsi = -1j*ham.dot(y) - hdecay;
		return Hpsi
	# -----------------------
	# make full hamiltonian:
	if loopover == "lambda0":
		lamb0 = o.lambin0[il];
		ham = o.ham1 + wv*lamb0*o.Hbsm + wv*lamb0**2*o.sft;
		#print('il, o.ld*lamb0 = ',il,o.ld,lamb0)
		psi0 = createpsi0(o.ld*lamb0);
		#print('n*wv*lamb**2 = ',n*wv*(o.ld*lamb0)**2)
	else:
		g = o.lambin0[il]/np.sqrt(o.n);
		# print('il,o.lambin0[il], np.sqrt(o.n), g = ',il,o.lambin0[il],np.sqrt(o.n),g)
		ham = o.ham1 + g*o.Hgsm;
		psi0 = createpsi0(o.ld*lambda0);
		# if o.n==1: print(ham.toarray())
	ham = ham.tocsr();

	t0 = 0;
	# start complex integrator:
	r = ode(fInteg).set_integrator('zvode', method='bdf')#,order=2)	r.set_initial_value(psi0, t0);
	# tlist = [0]; 
	corr = [1];# correlation fun at t0
	i = 1; # start from 1, t=0 already done!
	while r.successful() and i < ntmax: # r.t < tfddd: #
		tc = i*dt; # r.t+dt;# 
		print('xx = ',o.xx)
		psit = r.integrate(tc);
		# get the correlation function:
		corr.append(np.vdot(psi0,psit));
		# tlist.append(tc);
		# iterate the integrator; from current time 
		r.set_initial_value(psit, tc);
		i += 1;
	if not r.successful():
		print("getd: r.successful() = F !! stop?!");
	corr = np.array(corr); # tlist = np.array(tlist);
	del r
#--------------------------
	return corr # , tlist


#----------------------------------------------------------
def gettdRK4(il):
	# time evolves the state with only a photon present
	# psi0 = |P(t=0)> = a^dag|0x,0p,0v> = [1,0,0,......,0];
	#--------------------------		
	dth = 0.5*dt; dt6 = dt/6;
	kappa2 = kappa/2.0; gamma2 = gamma/2.0;
	def yprime(y):
		if dkapa:
			hdecay = np.concatenate((kappa2*y[range(n1)],gamma2*y[range(n1,ntot)]));
		else:
			hdecay = kappa2*y;
		Hpsi = -1j*ham.dot(y) - hdecay;
		return Hpsi
	#--------------------------
	def RK4Step(y):
		k1 = yprime(y         )
		k2 = yprime(y + k1*dth)
		ks = k1+2*k2; k1=0; # defining ks saves 2 arrays, 3 instead of 5.
		k3 = yprime(y + k2*dth)
		ks += 2*k3; k2=0;
		k4 = yprime(y + k3*dt )
		ks += k4; k3=0;k4=0;
		return y+ks*dt6; #y + (k1+2*k2+2*k3+k4)*dt/6.
 	#--------------------------
	# make full hamiltonian:
	if loopover == "lambda0":
		lamb0 = o.lambin0[il];
		ham = o.ham1 + wv*lamb0*o.Hbsm + wv*lamb0**2*o.sft;
		#print('o.sft.data = ',o.sft.data)
		if nlmax ==1:
			o.Hbsm = []; o.sft=[];
		psi0 = createpsi0(o.ld*lamb0);
	else:
		g = o.lambin0[il]/np.sqrt(o.n);
		ham = o.ham1 + g*o.Hgsm;
		if nlmax ==1:
			o.Hgsm = [];
		psi0 = createpsi0(o.ld*lambda0);
	ham = ham.tocsr();
	# integrate
	t0 = 0; corr = [1]; i = 1; # start from 1, t=0 already done!
	psit = psi0;
	while i < ntmax:
		if i%printstep == 0:
			print(' time evolution step = '+str(i)+'/'+str(ntmax))
		tc = i*dt; 
		psit = RK4Step(psit);
		corr.append(np.vdot(psi0,psit));
		i += 1;
	corr = np.array(corr);
	return corr
#----------------------------



def gettdfort(il):
	import tcorr # fortran routine
	# time evolves the state with only a photon present
	# psi0 = |P(t=0)> = a^dag|0x,0p,0v> = [1,0,0,......,0];
	#--------------------------		
	#dth = dt/2; dt6 = dt/6;
	kappa2 = kappa/2.0; gamma2 = gamma/2.0;
 	#--------------------------
	# make full hamiltonian:
	if loopover == "lambda0":
		lamb0 = o.lambin0[il];
		ham = o.ham1 + wv*lamb0*o.Hbsm + wv*lamb0**2*o.sft;
		#print('o.sft.data = ',o.sft.data)
		if nlmax ==1:
			o.Hbsm = []; o.sft=[];
		psi0 = createpsi0(o.ld*lamb0);
	else:
		g = o.lambin0[il]/np.sqrt(o.n);
		ham = o.ham1 + g*o.Hgsm;
		if nlmax ==1:
			o.Hgsm = [];
		psi0 = createpsi0(o.ld*lambda0);
	ham = ham.tocsr();
	# create arguments for fortran routine
	Val=ham.data;# CSR format data array of the matrix
	Col=ham.indices;# CSR format index array of the matrix
	RowPtr=ham.indptr;# CSR format index pointer array of the matrix
	nnz, = Val.shape
	nrp, = RowPtr.shape	
	# shift index for fortran:
	Col += 1; RowPtr += 1;
	# integrate
	corr = tcorr.tcorr(n1,ntot,nnz,nrp,printstep,Val,Col,RowPtr, psi0,dt,ntmax,kappa2,gamma2)
	return corr
#----------------------------





# our spectrum has only nonzero peaks in a small reagion
# FFT/IFFT not good because they spread the points on freq axis evenly and only a few ~ 50 points among ~10k lie within relevent region where we get absorption peaks, i.e., ~ 2.0 around wx.
# So, I simply get the FT in the desired freq range with desired point density for a smooth spectrum. 
# Possible issue:
#	Due to the correlation function not well decayed to zero:
#		we can get an abrupt change at the end, creating a fast oscillating feature in the FT. 
#		1. for small system sizes, we can evolve the system for longer times till the correlation function decays well enough to zero.
#		2. for large system sizes, we can just filter this unphysical feature by taking fft and taking ifft with only low freq components.
#----------------------------------------------------------
def getFourierTransform(Xt,ti,dd,wlist,arrow):
	Xw = []; ftnorm = 1/np.sqrt(2*pi);
	p = 1j;
	if arrow == 'B':
		p = -1j;
	for w in wlist:
		xw = 0; 	t= ti;
		for xt in Xt:
			xw += xt*exp(p*w*t);
			t += dd;
		Xw.append(xw*dd*ftnorm);
	Xw = np.array(Xw);
	return Xw
#--------------------------
#---------------------------------
def fwriteCorr(il,tlist, corr, fnametd):
	# save correlation vs time
	ntmax = len(tlist);
	corrOUT = np.zeros((ntmax,3));
	corrOUT[:,0] = tlist;
	corrOUT[:,1] = np.real(corr);
	corrOUT[:,2] = np.imag(corr);
	f=open(fnametd,'ab');
	if loopover =="lambda0":
		lamb0 = o.lambin0[il];
		wr0 = wr;
	else:
		wr0 = o.lambin0[il];
		lamb0 = lambda0;
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, corrOUT,fmt='%15.10f %15.10f %15.10f', delimiter=' ', header=header,comments='#')
	f.close();
	f=open(fnametd,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return None
#--------------------------
def fwriteFT(il,wlist, Gw, GR, fnametd, ntmax):
	# save FT of correlation
	absOUT = np.zeros((nwmax,4));
	absOUT[:,0] = wlist;
	absOUT[:,1] = np.real(Gw);
	# 2*kl*np.imag(absws)**2 + kappa*np.abs(absws)**2;
	absOUT[:,2] = 2*kl*np.real(Gw)**2 + kappa*np.abs(Gw)**2;
	absOUT[:,3] = -np.imag(GR); # Green from analytical
	f=open(fnametd,'ab');
	if loopover =="lambda0":
		lamb0 = o.lambin0[il];
		wr0 = wr;
	else:
		wr0 = o.lambin0[il];
		lamb0 = lambda0;
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr0)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, absOUT,fmt='%15.10f %15.10f %15.10f %15.10f', delimiter=' ', header=header,comments='#')
	f.close();
	f=open(fnametd,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return None

def fwriteFT2(il,wlist, Gw, GR, GR2, fnametd, ntmax):
	# save FT of correlation
	absOUT = np.zeros((nwmax,4));
	absOUT[:,0] = wlist;
	absOUT[:,1] = np.real(Gw);
	# 2*kl*np.imag(absws)**2 + kappa*np.abs(absws)**2;
	#absOUT[:,2] = 2*kl*np.real(Gw)**2 + kappa*np.abs(Gw)**2;
	absOUT[:,2] = -np.imag(GR2);# truncated then transformed
	absOUT[:,3] = -np.imag(GR); # Green from analytical (transformed then trucated)
	f=open(fnametd,'ab');
	if loopover =="lambda0":
		lamb0 = o.lambin0[il];
		wr0 = wr;
	else:
		wr0 = o.lambin0[il];
		lamb0 = lambda0;
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr0)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, absOUT,fmt='%15.10f %15.10f %15.10f %15.10f', delimiter=' ', header=header,comments='#')
	f.close();
	f=open(fnametd,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return None

#--------------------------------------------
# Analytical green function from JK notes
# -------------------------------------------
def getFrankCondonEtc(l0):
	res=[]; 	efac=[]; #mx = 2*o.mx
	if abs(l0) < 1e-8:
		for i in range(mx+1):
			if i==0: x=1;
			else: x=0;
			res.append(x);
			efac.append(i*wv + wx -1j*gamma/2);
		return res,efac
	cont=decimal.Context(prec=15, Emax=999, clamp=1);
	decl0 = decimal.Decimal.from_float(l0);
	for i in range(mx+1):
		#x = np.exp(-l0**2)*l0**(2*i)/factorial(i);
		n0 = decimal.Decimal.from_float(factorial(i));
		fc0 = cont.power(decl0,2*i);
		x = np.exp(-l0**2)*float(fc0/n0);
		res.append(x);
		efac.append(i*wv + wx - wv*l0**2 -1j*gamma/2); # -wv*l0**2 polaron transform term not included in JK notes
	return res,efac
# -------------------------------------------
def Green(wlist,l0,wr):
	#mx = 2*o.mx
	# ------------------
	# getFrankCondonEtc()
	FC, efac = getFrankCondonEtc(l0);# Frank-Condon factors etc
	# ------------------
	#print(wr,wx,wc,wv,l0)

	wr2 = wr*wr;
	GRlist = [];
	for w in wlist:
		SelfE = 0;
		for i in range(mx+1):
			SelfE -= wr2*FC[i]/(w - efac[i]);
		GR = 1/(w - wc + 1j*kappa/2 + SelfE);
		GRlist.append(GR);
	GRlist = np.array(GRlist);
	ftnorm = 1/np.sqrt(2*pi);
	return GRlist*ftnorm
# -------------------------------------------
def fcorrft(il):
	# calculate correlation
	memtime('corr');
	if usenumpy:
		corr = gettdRK4(il); # numpy
		memtime('end-numpy');
	else:
		corr = gettdfort(il);# fortran
		memtime('end-fort');
	E1, E2 = e1+wx, e2+wx;	
	#print('using: e1,e2 = ',e1,e2)
	wlist = np.linspace(E1, E2,nwmax);
	Gw = getFourierTransform(corr,0,dt,wlist,'F');
	# green function based results:
	# wlist = wlist - wx ; # shift freq axis
	if loopover =="lambda0":
		lamb0 = o.lambin0[il];
		wr = o.wr;
	else:
		lamb0 = lambda0;
		wr = o.lambin0[il];
	GR = Green(wlist,lamb0,wr); 	# analytical Green function
	GR2 = Green2(wlist,lamb0,wr); 	# analytical Green function

	return Gw, GR, GR2, corr
# -------------------------------------------
# create psi0: important case is when ld>0, coherent state for every site
def createpsi0(lam0):
	# psi0 already stored? when not stores: o.psi0=[];
	if len(o.psi0)>0:
		return o.psi0	
	# create:
	if ld>0 and lam0 > 0: # for displaced basis
		psi0 = np.zeros(ntot);
		psi0p=[]; row=[]; col=[];
		cont=decimal.Context(prec=15, Emax=999, clamp=1);
		decl0 = decimal.Decimal.from_float(lam0);
		for i in range(n1fsym):
			ii= int(o.Fact1l[i]);
			Nv = int(o.Nv1l[i]);
			xa = cont.power(decl0,Nv);
			y = decimal.Decimal(ii).sqrt();
			z = xa/y
			x = float(z)* exp(-o.n*lam0**2/2) * np.sqrt(o.Norm1l[i]);
			if x> 1e-12: # remove very small numbers to avoid numpy/scipy crash
				psi0[i] = x;
		# ----------- 'extra' basis --------------
		mn = n1fsym - (m+1)*n2;
		for jj in range(n2):
			ii= int(o.Fact1l[jj]);
			Nv = int(o.Nv1l[jj]); 
			xa = cont.power(decl0,Nv);
			y = decimal.Decimal(ii).sqrt();
			z = xa/y;
			x = float(z)* exp(-o.n*lam0**2/2) * np.sqrt(o.Norm2l[jj]);
			for mj in range(m+1,mx+1):
				factmj = int(factorial(mj));
				factmjsq = decimal.Decimal(factmj).sqrt();
				xa = cont.power(decl0,mj);
				x *= float(xa/factmjsq);
				if x> 1e-12: # remove very small numbers to avoid numpy/scipy crash
					i = mn + mj*n2 + jj;
					psi0[i] = x;
		norm=np.linalg.norm(psi0)
		if norm > 0:
			psi0 = psi0/norm;
	else: # for undisplaced basis
		psi0 = np.zeros(ntot);
		psi0[0] = 1;
	#psi0 = psi0.tocsr();
	# store this to global array for reuse if wr loop
	if loopover=='wr':
		o.psi0 = psi0;
	return psi0

# -----------------------------------------------------
# absoprtion option 1: time evolution
# get photon correlation function, Fourier transform it. 
# -----------------------------------------------------
#if (absorption=='true' and td == 'true'):

def corrft():
	print(' ====> calculating absorption: time evolution ... ');
	# indices for projectors: can we work with sparse?
	# o.indPp = np.arange(0,n1); o.indPx = np.arange(n1,ntot);
	# ------------------------------
	# set o.ham1
	if loopover == 'lambda0':
		g= wr/np.sqrt(n);
		o.ham1 =wc*o.Hcsm + wx*o.Hxsm + g*o.Hgsm + wv*o.Hvsm
		o.Hcsm=[]; o.Hxsm=[]; o.Hgsm=[]; o.Hvsm=[]
		# wv*o.Hbsm + wv*lamb0**2*o.sft
		# ham.tocsr();# csr for fast matrix vector products.
	elif loopover == 'wr':
		# print(wc,wv,wv,lambda0)
		o.ham1 = wc*o.Hcsm + wx*o.Hxsm +wv*o.Hvsm +lambda0*wv*o.Hbsm + wv*lambda0**2*o.sft;
		o.Hcsm=[]; o.Hxsm=[]; o.Hbsm=[]; o.Hvsm=[]; o.sft= [];
	# ------------------------------
	# print(o.ham1.toarray())
	# Number of processes
	if nlmax<Np:
		Npl = nlmax;
	else:
		Npl = Np
	# print(nlmax,Np,Npl)
	# ------------------------------
	if nlmax == 1: results = [fcorrft(0)];
	else:
		# fresh pool
		pool=Pool(Npl);
		results = pool.map(fcorrft, range(nlmax));
		o.ham1 = []; o.Hbsm=[];o.sft=[]; o.Hgsm=[];
		pool.close();
	# ------------------------------
	#print('writing: e1,e2 = ',e1,e2)
	wlist = np.linspace(e1,e2,nwmax);
	tlist = np.linspace(0, tf, ntmax);
	print(" ntmax = ",ntmax);
	print(" nwft = ",nwmax);
	# ------------------------------
	il = 0;
	for Gw, GR, GR2, corr in results:
		# 1/2.5 makes GR peaks equal to the exact for lam=0 case.
		fwriteFT2(il, wlist, Gw, GR, GR2, dumy+'/abs-vs-w-td.txt', ntmax); # write FT file
		fwriteCorr(il, tlist, corr, dumy+'/corr-vs-t-td.txt');# write Correlation file
		il += 1;
	return
#------------------------------------------------------

def Green2(wlist,l0,wr):
	#mx = 2*o.mx
	h = np.zeros((mx+1,mx+1));
	for i in range(mx+1):
		h[i,i] = i;	
	for i in range(mx):
		h[i,i+1] = np.sqrt(i+1)*l0;
	for i in range(1,mx+1):
		h[i,i-1] = np.sqrt(i)*l0;

	es,evs = linalg.eig(h); del h;
	#print(evs)
	#print('-----')
	evs2 = evs[0,:]**2;
	#print(evs2)
	
	#print(linalg.norm(evs,axis=0))

	#print(wr,wx,wc,wv,l0)
	wr2 = wr*wr;
	wxc = wx - 1j*gamma/2; wcc = wc - 1j*kappa/2;
	GRlist = [];
	for w in wlist:
		SelfE = 0;
		for i in range(mx+1):
			SelfE -= wr2*evs2[i]/(w - wxc - wv*es[i]);
		GR = 1/(w - wcc + SelfE);
		GRlist.append(GR);
	GRlist = np.array(GRlist);
	ftnorm = 1/np.sqrt(2*pi);
	return GRlist*ftnorm








