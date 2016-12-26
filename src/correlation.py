

import globalvariables as o
import numpy as np;
from math import factorial
from scipy.integrate import ode
from numpy import pi, exp
from scipy.sparse import coo_matrix, diags
import matplotlib.pyplot as plt
from multiprocessing import Pool
import decimal

# define local to avoid writing o. every time!
# at risk of setting these at the time of first import of this module only.
##### any update in the values of these variables will not be considered by the functions of this module, because they use the first-import-time's set values.
########### consider this if loops over lambda/wr/n for absorption calculations is tried.
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2,ntot = o.n1,o.n2,o.ntot;
n1fsym = o.n1fsym;
ld = o.ld;
tolr, itermax = o.tolr, o.itermax 


gamma, kappa =	 o.gamma, o.kappa
if abs(gamma-kappa) > 1e-5: dkapa = 1;
else: dkapa = 0;

e1,e2,dt,tf = 	o.e1, o.e2, o.dt, o.tf
nwmax,ntmax =	o.nwmax, o.ntmax
mg = o.mg
show = 	o.show

kl= kappa/2;

lmin,lmax, nlmax = o.lmin, o.lmax,  o.nlmax;
loopover= o.loopover;
lambda0= o.lambda0;


# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

#----------------------------------
#	absorption using td: functions
	#	steps:
		# start with |P(t=0)> = a^dag|0x,0p,0v>, 
		# evolve to get |P(t)> & get the correlation function <P(0)|P(t)>
		# Fourier transform to get the absorption
#----------------------------------
def getwlist(ntmax,dt):
	if ntmax%2 == 0:
		#f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)
		f = np.linspace(-ntmax/2,ntmax/2-1,ntmax);
	else:
		#f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
		f = np.linspace(-(ntmax-1)/2,(ntmax-1)/2,ntmax);
	f *= 1/(dt*ntmax);
	return f;
#--------------------------
def gettd(il):
	# time evolves the state with only a photon present
	# psi0 = |P(t=0)> = a^dag|0x,0p,0v> = [1,0,0,......,0];
	#--------------------------
	# making internal function
	# to avoid making list of H for loopover lam/wr
	#--------------------------
	def fInteg(t, y):
		# gives RHS of td-schrodinger equation with losses:
		# gives [-iota*H - (kappa*PhotProjector + gamma*ExcProjector)]psi(t)
		if dkapa: hdecay = np.concatenate( (kappa*y[range(n1)],gamma*y[range(n1,ntot)]) );
		else: hdecay = kappa*y	
		Hpsi = -1j*ham.dot(y) - hdecay;
		return Hpsi
	# -----------------------
	# make full hamiltonian:
	if loopover == "lambda0":
		lamb0 = o.lambin0[il];
		ham = o.ham1 + wv*lamb0*o.Hbsm + wv*lamb0**2*o.sft;
		psi0 = createpsi0(o.ld*lamb0);
		print('n*wv*lamb**2 = ',n*wv*(o.ld*lamb0)**2)
	else:
		g = o.lambin0[il]/np.sqrt(o.n);
		# print('il,o.lambin0[il], np.sqrt(o.n), g = ',il,o.lambin0[il],np.sqrt(o.n),g)
		ham = o.ham1 + g*o.Hgsm;
		psi0 = createpsi0(o.ld*lambda0);
		# if o.n==1: print(ham.toarray())
	ham = ham.tocsr();

	# psi0 =  Psi0fromPhotgs(ham); # use only if 
			# issues with large factorial handling etc
			#  are not resolved with decimal module
	if 0:
		# print(' handmade: Psi0 = --------------------')
		# print(psi0)
		nplt = 400
		# print('hand: ',psi0[51:55])
		# print('ind51: Nv,Pi =',o.Nv1l[51],o.Norm1l[51:55])
		plt.plot(range(nplt),psi0[0:nplt],'sr',ms=10)
		psi0 =  Psi0fromPhotgs(ham);
		# print('calc: ',psi0[51:55])
		print('Psi0gs = --------------------')
		plt.plot(range(nplt),np.abs(psi0[0:nplt]),'og',ms=5)
		plt.show()
		_exit()

	t0 = 0;
	# start complex integrator:
	r = ode(fInteg).set_integrator('zvode', method='bdf')
	r.set_initial_value(psi0, t0);
	# tlist = [0]; 
	corr = [1];# correlation fun at t0
	i = 1; # start from 1, t=0 already done!
	while r.successful() and i < ntmax: # r.t < tfddd: #
		tc = i*dt; # r.t+dt;# 
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
#--------------------------
def getFT():
	# Take Fourier Transform: n= 2*ntmax
	ntime= 1*ntmax + 0;
	Gw = np.fft.ifft(corr,n=ntime, norm="ortho");
	Gw = np.fft.fftshift(Gw); # for correct order of freq
	#Gw = fft(corr); # scipy
	#wlist = np.fft.fftfreq(ntmax,d=dt);
	wlist = getwlist(ntmax,dt); # correct ordered freq
	wlist *= 2*pi; # to get right scale
	# do we need to calc this full wlist fist?
	# select the relevant freq window
	i1,i2 = getFrewWindowCoor(wlist,e1+wx,e2+wx); #+wx because, unshifted scale
	wlist = wlist[i1:i2];
	Gw = Gw[i1:i2]	;
	return Gw, wlist
#---------------------------------
def getFrewWindowCoor(wlist,e1,e2):
	# check if e1,e2 window is bigger?
	if e1 < wlist[0] or e2 > wlist[ntmax-1]:
		e1 = wlist[0]; i1=0;
		e2 = wlist[ntmax-1]; i2=ntmax;
		return i1,i2
	# start search from centre of array:
	ic = ntmax//2;
	i1 = 0; i2 = 0
	while wlist[ic-i1] > e1:
		i1 += 1;
	i1 = ic - i1;
	while wlist[ic+i2] < e2:
		i2 += 1;
	i2 = ic + i2;
	return i1,i2
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
#--------------------------------------------
# Analytical green function from JK notes
# -------------------------------------------
def getFrankCondonEtc(l0):
	res=[]; 	efac=[];
	for i in range(mg+1):
		x = np.exp(-l0**2)*l0**(2*i)/factorial(i);
		res.append(x)
		efac.append(i*wv -wv*l0**2 -1j*gamma/2); # -wv*l0**2 polaron transform term not included in JK notes
	return res,efac
# ---------------------
def Green(wlist,l0):
	# ------------------
	# getFrankCondonEtc()
	FC, efac = getFrankCondonEtc(l0);# Frank-Condon factors etc
	# ------------------
	wr2 = wr*wr; delta = wc-wx;
	GRlist = [];
	for w in wlist:
		SelfE = 0;
		for i in range(mg+1):
			SelfE -= wr2*FC[i]/(w - efac[i]);
		GR = 1/(w +1j*kappa/2 - delta + SelfE);
		GRlist.append(GR);
	GRlist = np.array(GRlist);
	return GRlist
# -------------------------------------------

def fcorrft(il):
	corr = gettd(il); # calculate correlation
	E1, E2 = e1+wx, e2+wx;		
	wlist = np.linspace(E1, E2,nwmax);
	Gw = getFourierTransform(corr,0,dt,wlist,'F');
	# green function based results:
	# wlist = wlist - wx ; # shift freq axis
	if loopover =="lambda0": lamb0 = o.lambin0[il];
	else: lamb0 = lambda0;
	GR = Green(wlist-wx,lamb0); 	# analytical Green function
	return Gw, GR, corr


# -----------------------------------------------------
# calculate the initial state, no vibrations, just a photon.
# in displaced basis, this is a coherent state for all sites
# we can construct it using createpsi0() below
# or calculate it using Psi0fromPhotgs() as the gs of 
#  photon block with zeros padded at the end to make dim = ntot.
# ----------------------------------
# These two should match only for a large vib cutoff because
#	the constructed one is for infinite cutoff whereas
# 	the calculated one is obviously for a finite cutoff.
# -----------------------------------------------------

# -----------------------------------------------------
def Psi0fromPhotgs(Hfull):
	# psi0 already stored? when not stores: o.psi0=[];
	if len(o.psi0)>0: 	
		return o.psi0	
	# create:
	if ld>0: # for displaced basis
		Hphot=Hfull[0:n1,0:n1];
		from scipy.sparse.linalg import lobpcg
		ev0 = np.random.rand(n1,1);
		evalu, evec = lobpcg(A=Hphot, X=ev0, tol=10e-10, maxiter=150,largest=False, verbosityLevel=0);	
		print('Hphot: Eval = ',evalu);
		psi0 = np.zeros(ntot);
		psi0[0:n1] = evec; # evec.T[0];
	else: # for undisplaced basis
		psi0 = np.zeros(ntot);
		psi0[0] = 1;
	# store this to global array for reuse if wr loop
	if loopover=='wr':
		o.psi0 = psi0;
	return psi0
#--------------------------------------


	

# create psi0: important case is when ld>0, coherent state for every site
def createpsi0(lam0):
	# psi0 already stored? when not stores: o.psi0=[];
	if len(o.psi0)>0:
		return o.psi0	
	# create:
	if ld>0: # for displaced basis
		psi0 = np.zeros(ntot);
		psi0p=[]; row=[]; col=[];
		cont=decimal.Context(prec=15, Emax=999, clamp=1);
		decl0 = decimal.Decimal.from_float(lam0);
		for i in range(n1fsym):
			ii= int(o.Fact1l[i]);
			Nv = int(o.Nv1l[i]);
			xa = cont.power(decl0,Nv);
			#print('------- xa = ',round(xa,5))
			if ii <0: 
				print('---- i ,ii fit kuj??? = ',i,ii)
			y = decimal.Decimal(ii).sqrt();
			z = xa/y
			# print('-----------Decimal(50)-------')
			# print('x,y,z = ',x,y,z)	
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
		# ----------- 'extra' basis --------------
	else: # for undisplaced basis
		psi0 = np.zeros(ntot);
		psi0[1] = 1;
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
	# loop parameter: lambbda_0 or wr
	lambin0 = np.linspace(lmin,lmax, nlmax);
	o.lambin0 = lambin0;
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
	
	wlist = np.linspace(e1,e2,nwmax);
	tlist = np.linspace(0, tf, ntmax);
	print(" ntmax = ",ntmax);
	print(" nwft = ",nwmax);
	# ------------------------------
	il = 0;
	for Gw, GR, corr in results:
		#print(np.abs(Gw))
		fwriteFT(il, wlist, Gw, GR, dumy+'/abs-vs-w-td.txt', ntmax); # write FT file
		fwriteCorr(il, tlist, corr, dumy+'/corr-vs-t-td.txt');# write Correlation file
		il += 1;
	# ------------------------------
	if 0: #show:
	
		#print(GR.shape)
		plt.plot(wlist, np.real(Gw),'-r',lw=1,label='$Re[G(w)]$')
		#plt.plot(wlist, np.real(Gwtw),'-g',lw=1,label='Gwtw')
		plt.plot(wlist, -np.imag(GR),'-g',lw=1,label='$Im[GR(w)]$')
		plt.legend()
		plt.figure()
		plt.plot(tlist, corr.real,'-r',lw=1,label='$Re[G(t)]$');
		#plt.plot(tlist, Gwt.real,'-b',label='Gwt');
		plt.plot(tlist, corr.imag,'g',lw=1,label='$Im[G(t)]$')
		plt.legend()
		plt.show()
	return
#------------------------------------------------------











