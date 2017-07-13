from qutip import *
from numpy import linspace, meshgrid, sqrt
from scipy.sparse.linalg import lobpcg
from scipy.sparse import identity
import numpy as np

import os,sys
os.makedirs('fullHilbert/data/n-all', exist_ok=True)

# path to input parameter file
path_param = os.getcwd();
# add this path to python search path
sys.path.insert(0, path_param)


#from parameters import *
import param
Om = param.wv;
omc = param.wc;
om0 = param.wx;
if param.loopover == 'lambda0':
	omR = param.wr;
	lmin = param.lmin
	lmax = param.lmax
	nlmax = param.nlmax;
	lamlist = np.linspace(lmin,lmax,nlmax);
else:
	print(' use loopover=lambda0 please.... ')
M = param.mlist[0][1]-1;
print('param.mlist = ',param.mlist)
print('m,wr,wv,wc,wx = ',M-1, omR,Om,omc,om0)
print('lamlist = ',lamlist)

ntot = 3*M**2;
n1 = M**2;
if ntot > 30:
	nstates=ntot-5;
else:
	nstates=ntot-2;

print('ntot, nstates = ',ntot,' ',nstates)

def mkH(lam0):
	b = destroy(M)
	Ib = qeye(M)
	Is = qeye(3)
	phot0 =basis(3,0)
	phot = phot0*phot0.dag()
	S10 = basis(3,1)
	S1 = S10*S10.dag()
	S20 = basis(3,2)
	S2 = S20*S20.dag()
	H = omc*tensor(phot, Ib, Ib) +om0*tensor(S1+S2, Ib, Ib) 
	H += omR/sqrt(2)*tensor(phot0*S10.dag() + phot0*S20.dag() + S10*phot0.dag() + S20*phot0.dag(), Ib, Ib)
	H += Om*tensor(Is, b.dag()*b,Ib) + Om*tensor(Is, Ib, b.dag()*b) + Om*lam0*tensor(S1, b.dag()+b, Ib) +  Om*lam0*tensor(S2, Ib, b.dag()+b)
	#print('H.shape = ', H.shape)
	return H.data

def elb(l0):
	ret = 1/2*(omc+om0-Om*l0**2) - np.sqrt(((om0-omc-Om*l0**2)/2)**2+omR**2)
	return ret

def diagonlise(H):
	iden = identity(ntot).toarray();
	ev0 = np.random.rand(ntot, nstates) ; # giving 1 vec means k=1, lowest eigenpair.
	enlb=elb(lam0);
	#print(H.shape, iden.shape)
	ham = H - enlb*iden;
	evalu, evec = lobpcg(A=ham, X=ev0, tol=10e-8, maxiter=200,largest=False,verbosityLevel=1);	
	elp = evalu+enlb;
	return elp-om0,evec
	
def writeout(lam0,evs,mep,mex,cos2th):
	# combine the results for writing output file:
	absOUT = np.zeros((nstates,4));
	absOUT[:,0] = evs;
	absOUT[:,1] = mep;
	absOUT[:,2] = mex;
	absOUT[:,3] = cos2th;
	fabsorption = "fullHilbert/data/n-all/absorption.txt"
	f=open(fabsorption,'ab');
	header = " "+str(nstates)+" "+str(2)+" "+str(M-1)+" "+str(M-1);
	header += " "+str(omR)+" "+str(lam0)+" "+str(omc)+" "+str(om0)+" "+str(Om);
	np.savetxt(f, absOUT,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ',	header=header,comments='#')
	f.close();
	f=open(fabsorption,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return


evalues = []; matelem = []; matelemx=[];
for lam0 in lamlist:
	print('lam0 = '+str(lam0)+' .... ')
	H =mkH(lam0);
	evalu, evec = diagonlise(H);
	evs = evalu;
	mep = np.abs(evec[0,:])**2;
	mex = np.abs(evec[n1,:])**2;
	cos2th = np.sum(np.abs(evec[0:n1,:])**2,axis=0)
	writeout(lam0,evs,mep,mex,cos2th)

	evalues.append(evalu);
	matelem.append(mep);
	

evalues = np.array(evalues);
evalues.dump('es-full');

matelem = np.array(matelem);
matelem.dump('matelem-full');

print(' completed ..... :) ')


	
	

