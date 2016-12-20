

import globalvariables as o
import numpy as np;
from math import factorial
from scipy.integrate import ode
from numpy import pi, exp
from scipy.sparse import coo_matrix, diags
import matplotlib.pyplot as plt


# define local to avoid writing o. every time!
# at risk of setting these at the time of first import of this module only.
##### any update in the values of these variables will not be considered by the functions of this module, because they use the first-import-time's set values.
########### consider this if loops over lambda/wr/n for absorption calculations is tried.
n,m,mx = o.n, o.m, o.mx
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2,n3,ntot = o.n1,o.n2,o.n3,o.ntot;
lamb0 = o.lamb0;

gamma, kappa =	 o.gamma, o.kappa
e1,e2,dt,tf = 	o.e1, o.e2, o.dt, o.tf
nwmax =	o.nwmax
mg = o.mg
show = 	o.show

kl= kappa/2;

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
def fInteg(t, y):
	# gives RHS of td-schrodinger equation with losses:
	# uses global: ham1, kappa, gamma
	# takes: psi(t);
	# gives [-iota*H - (kappa*PhotProjector + gamma*ExcProjector)]psi(t)
	hdecay = np.concatenate( (kappa*y[range(n1)],gamma*y[range(n1,ntot)]) );
	Hpsi = -1j*o.ham1.dot(y) - hdecay;
	return Hpsi
#--------------------------
def gettd():
	# time evolves the state with only a photon present
	# psi0 = |P(t=0)> = a^dag|0x,0p,0v> = [1,0,0,......,0];
	row=[0];col=[0];dat=[1];
	psi0=coo_matrix((dat, (row, col)), shape=(1,ntot));
	# psi0 = psi0.tocsr();
	psi0 = psi0.toarray()[0];
	t0 = 0;
	# start complex integrator:
	r = ode(fInteg).set_integrator('zvode', method='bdf')
	r.set_initial_value(psi0, t0);
	tlist = [0]; corr = [1];# correlation fun at t0
	while r.successful() and r.t < tf:
		tc = r.t+dt; 
		psit = r.integrate(tc);
		# get the correlation function:
		corr.append(np.vdot(psi0,psit));
		tlist.append(tc);
		# iterate the integrator; from current time 
		r.set_initial_value(psit, tc);
	corr = np.array(corr); tlist = np.array(tlist);
	return corr, tlist
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
def fwriteCorr(tlist, corr, fnametd):
	# save correlation vs time
	ntmax = len(tlist);
	corrOUT = np.zeros((ntmax,3));
	corrOUT[:,0] = tlist;
	corrOUT[:,1] = np.real(corr);
	corrOUT[:,2] = np.imag(corr);
	f=open(fnametd,'ab');
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, corrOUT,fmt='%15.10f%15.10f%15.10f', delimiter=' ', header=header,comments='#')
	f.close();
	f=open(fnametd,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return None
#--------------------------
def fwriteFT(wlist, Gw, GR, fnametd, ntmax):
	# save FT of correlation
	absOUT = np.zeros((nwmax,4));
	absOUT[:,0] = wlist;
	absOUT[:,1] = np.real(Gw);
	# 2*kl*np.imag(absws)**2 + kappa*np.abs(absws)**2;
	absOUT[:,2] = 2*kl*np.real(Gw)**2 + kappa*np.abs(Gw)**2;
	absOUT[:,3] = -np.imag(GR); # Green from analytical
	f=open(fnametd,'ab');
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, absOUT,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ', header=header,comments='#')
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















# -----------------------------------------------------
# absoprtion option 1: time evolution
# get photon correlation function, Fourier transform it. 
# -----------------------------------------------------
#if (absorption=='true' and td == 'true'):
def corrft():
	print(' ====> calculating absorption: time evolution ... ');
	# indices for projectors: can we work with sparse?
	# o.indPp = np.arange(0,n1); o.indPx = np.arange(n1,ntot);
	
	# Hamiltonian:
	g= np.sqrt(n); g= wr/g;

	ham = wc*o.Hcsm + wx*o.Hxsm + g*o.Hgsm + wv*o.Hvsm;
	ham += lamb0*wv*o.Hbsm + wv*lamb0**2*o.sft;
	
	o.ham1 = ham.tocsr(); 	# csr for fast matrix vector products.
	del o.Hcsm; del o.Hxsm; del o.Hgsm; del o.Hvsm;
	del o.Hbsm; del o.sft;

	corr, tlist = gettd(); # calculate correlation
	ntmax=len(tlist); # no. of time points
	print(" ntmax = ",ntmax);
	if 1:
		E1, E2 = e1+wx, e2+wx;		
		wlist = np.linspace(E1, E2,nwmax);
		if 0:
			cctf = corr[-1];
			corrDecay = []; tlistadd = [];
			tff= ntmax*dt;
			nsmooth = 100;
			for i in range(1,nsmooth):
				corrDecay.append(cctf*np.exp(-i));
				tff += dt
				tlistadd.append(tff);
			corrDecay = np.array(corrDecay);
			tlistadd = np.array(tlistadd);
			corr = np.concatenate((corr,corrDecay));
			tlist = np.concatenate((tlist,tlistadd));
			ntmax = ntmax+nsmooth-1;
		Gw = getFourierTransform(corr,0,dt,wlist,'F');
		if 0:
			# remove fast unphysical oscillations		
			dw = wlist[1] - wlist[0];
			Gwt = getFourierTransform(Gw,E1,dw,tlist,'B');
			Gwtw = getFourierTransform(Gwt,0,dt,wlist,'F');
	# Gw, wlist = getFT(); # get Fourier Transform of corr
	# Gw, wlist = getMyFT(corr); # get Fourier Transform of corr
	# green function based results:
	wlist = wlist - wx ; # shift freq axis
	GR = Green(wlist,lamb0); 	# analytical Green function
	fwriteFT(wlist, Gw, GR, dumy+'/abs-vs-w-td.txt', len(tlist)); # write FT file
	fwriteCorr(tlist, corr, dumy+'/corr-vs-t-td.txt');# write Correlation file

	print(" nwft = ",nwmax);
	
	if show:
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











