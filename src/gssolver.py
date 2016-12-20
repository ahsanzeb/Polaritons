
import globalvariables as o
import numpy as np
from scipy.sparse import identity
from multiprocessing import Pool
from eigsolver import fdiagl, fdiagw

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2,ntot = o.n1,o.n2,o.ntot;
detuning, lamb0, eshft, nstates = o.detuning, o.lamb0, o.eshft, o.nstates

lmin,lmax, nlmax = o.lmin, o.lmax,  o.nlmax;
loopover= o.loopover;
lambda0= o.lambda0;

#************************************************************
# ground state calculations with loop over lambda or wr
#************************************************************
def gssolve():
	print(' ====> calculating groundstate ... ')
	iden = identity(ntot);
	o.iden = iden.tocsc();
	o.ev0 = np.random.rand(ntot, nstates) ; # giving 1 vec means k=1, lowest eigenpair.

	lambin0 = np.linspace(lmin,lmax, nlmax);
	o.lambin0 = lambin0;

	lambin = []; il=-1;
	for lamb00 in lambin0:
		il=il+1
		lambin.append([il,lamb00])
	o.lambin = lambin;
	# Number of processes for diagonalisation loop
	if nlmax<Np:
		Npdiag = nlmax;
	else:
		Npdiag = Np

	if loopover == 'lambda0':
		print('      calculating eigenvector/values... loop over lambda0 ')
		g= np.sqrt(n); g= wr/g;
		if (detuning==1):
			print("      detuning != 0 ")
			print("      ham1 =wc*Hcsm + wx*Hxsm + g*Hgsm + wv*Hvsm")
			o.ham1 =wc*o.Hcsm + wx*o.Hxsm + g*o.Hgsm + wv*o.Hvsm
			del o.Hcsm # free memory
			del o.Hxsm
		else:
			print("      detuning = 0 ")
			print("      ham1 = g*Hgsm + wv*Hvsm")
			o.ham1 = g*o.Hgsm + wv*o.Hvsm

		del o.Hgsm
		del o.Hvsm
		pool=Pool(Npdiag) # this would know ham1 and other variables
		o.eigvv = pool.map(fdiagl, lambin);
		#eigvv = pool.map(fdiagold, lambin)
		# ssss=gc.collect()
		# print('main: after loop: gc.collect() = ',ssss)
		del o.ham1 # free memory
		del o.Hbsm
	else: # if param.loopover == 'wr'
		print('      calculating eigenvector/values... loop over wr ')
		if (detuning == 1 ):
			o.ham1 = wc*o.Hcsm + wx*o.Hxsm + wv*o.Hvsm + lambda0*wv*o.Hbsm + wv*lambda0**2*o.sft;
			del o.Hcsm # free memory
			del o.Hxsm
		else:
			o.ham1 = wv*o.Hvsm + lambda0*wv*o.Hbsm + wv*lambda0**2*o.sft;
		del o.Hbsm
		del o.Hvsm
		pool=Pool(Npdiag) # this would know ham1 and other variables
		o.eigvv = pool.map(fdiagw, lambin)
		del o.ham1 # free memory
		del o.Hgsm
	
	pool.close() # this pool does not know eigvv
	return



