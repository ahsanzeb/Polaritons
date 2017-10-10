
import globalvariables as o
import numpy as np
from scipy.sparse.linalg import lobpcg
import sys


# define local to avoid writing o. every time!
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
eshft, nstates = o.eshft, o.nstates
tolr, itermax = o.tolr, o.itermax 
justenergy = o.justenergy
photonfraction = o.photonfraction;
lambda0= o.lambda0;

# LOWER BOUND ON ENERGY: lambda
def elb(l0):
	ret = 1/2*(wc+wx-wv*l0**2) - np.sqrt(((wx-wc-wv*l0**2)/2)**2+wr**2)
	return ret

# LOWER BOUND ON ENERGY: wr
def elbw(wr):
	ret = 1/2*(wc+wx-wv*lambda0**2) - np.sqrt(((wx-wc-wv*lambda0**2)/2)**2+wr**2)
	return ret


def fdiagl(lamb):
	il=lamb[0]; lamb0= lamb[1];
	# with sigma=E_LowerBound, A = ham-sigma*I
	# becomes positive definite (ham symmetric --> A symmetric) 
	enlb=elb(lamb0) # lower bound on E_LP
	ham = o.ham1 + lamb0*wv*o.Hbsm - enlb*o.iden + wv*lamb0**2*o.sft;
	if 1:
		evalu, evec = lobpcg(A=ham, X=o.ev0, tol=tolr, maxiter=itermax,largest=False,verbosityLevel=0);	
	#evalu,evec = scipy.sparse.linalg.eigsh(ham, k=nstates, which='SA',return_eigenvectors=True)
		# x = evec*evec;
		# print("norms:",np.sum(x,0))
		# print("norms:",np.sum(x,1))
	elif 1:
		evalu, evec = lobpcg(A=ham, X=o.ev0, tol=tolr, maxiter=itermax,largest=False,verbosityLevel=0)
	elif 0:
		evalu,evec = scipy.sparse.linalg.eigsh(ham, k=nstates, which='SA',return_eigenvectors=True)
		print('-----------')
		print(evalu)
	elif 1:
		evalu, evec, iiter,resNorm = lobpcg(A=ham, X=o.ev0, tol=tolr, maxiter=itermax,largest=False,verbosityLevel=1)
		if (resNorm > 10*tolr):
			print("il = ", il)
			print("resNorm > 10*tolr: resNorm,tolr= ",resNorm,tolr)	
			print("maxiter, iterused = ",itermax,iiter)
			print(" eigenvalues, vectors not converged!")

	elp = evalu+enlb + eshft;

	if nstates <= 10:
		print('il: evalu, enlb, eshft,elp=',il,evalu, enlb, eshft,elp)
		#print("il, lam0, Elp = ",il,lamb0, elp);
		sys.stdout.flush()
	if (justenergy and not photonfraction):
		return il, elp
	else:
		return il, elp, evec
#	gc.collect()
#****************************************************************
# diagonalisation: for wr loop. sigma=?
def fdiagw(lamb):
	il=lamb[0]; wr = lamb[1];
	g= wr/np.sqrt(o.n);
	# with sigma=E_LowerBound, A = ham-sigma*I
	# becomes positive definite (ham symmetric --> A symmetric) 
	enlb=elbw(wr) # lower bound on E_LP
	ham = o.ham1 + g*o.Hgsm - enlb*o.iden
	itermax = 100; tolr=1e-8; 
	evalu, evec = lobpcg(A=ham, X=o.ev0, tol=tolr, maxiter=itermax,largest=False,verbosityLevel=0)
	elp = evalu+enlb + eshft;
	if nstates <= 5:
		print("il, wr, Elp = ",il,wr, elp)
		sys.stdout.flush()	
	if (justenergy and not photonfraction):
		return il, elp
	else:
		return il, elp, evec
#****************************************************************
