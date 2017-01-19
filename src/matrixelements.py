


import globalvariables as o
from eigsolver import fdiagl, fdiagw
import numpy as np
from scipy.sparse import identity
from multiprocessing import Pool


# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2,ntot = o.n1,o.n2, o.ntot;
detuning, eshft, nstates = o.detuning, o.eshft, o.nstates

nlmax = o.nlmax;
loopover= o.loopover;
lambda0= o.lambda0;
lambin0 = o.lambin0;

# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

def getCos2th(i):
	cos2th = 0;
	# first n1 basis have a photon
	for p in range(n1): 
		cos2th += abs(o.evecs[p,i])**2
	return cos2th
#--------------------------
def phot_fraction(i):
	cos2th = 0;
	#evec = o.evecs[i];
	# first n1 basis have a photon
	for p in range(n1): 
		cos2th += abs(o.evecs[i][p])**2
	return cos2th
#--------------------------

# Herrera & Spano Absorption:
# Dipole matrix elements between eigenstates of HTC model in 1 excitation space (the one we calculate in the code).
def EiEjmatelem(i):
		matelem = np.zeros((nstates));
		eveci = o.evecs[n1:ntot,i]; # i-th eigenstate's exciton part
		for jj in range(n2):
			for mj in range(0,m+1):
				# n1+ mj*n2 + jj translated by -n1; exciton block basis index
				j= mj*n2 + jj; 
				ii = o.map21[jj,mj]; # photon block basis index
				matelem += o.evecs[ii,:]*eveci[j];
			# extra basis: diag in mj,jj
			for mj in range(m+1,mx+1): 
				# coor in exciton block
				j = mj*n2 + jj;
				# coor in photon block
				ii = n1fsym +(mj-m-1)*n2 + jj;
				matelem += o.evecs[ii,:]*eveci[j];
		return np.sum(np.abs(matelem)**2)
#--------------------------


# absoprtion option 2: get eigenstates and find matrix elements
# 	calc eigpairs in 1 excitation subspace and
# 	photon absorption matrix elements from 0 excitation subspace
# 	store for postprocessing, use util prog to 
# 	get absorption for whatever decay rates
# -----------------------------------------------------
# Calculate absorption spectrum related quantities:
# -----------------------------------------------------
# if (absorption=='true' and td != 'true'):


def fmatelem():
	print(' ====> calculating matrix elements ... ')
	# print('       calculating absorption: frequency space ... ');
	iden = identity(ntot);
	o.iden = iden.tocsc();
	o.ev0 = np.random.rand(ntot, nstates) ; # giving 1 vec means k=1, lowest eigenpair.

	# Hamiltonian:
	if loopover == 'lambda0':
		fsolve = fdiagl; # for diagonalisation
		g= wr/np.sqrt(n);
		if (detuning==1):
			o.ham1 =wc*o.Hcsm + wx*o.Hxsm + g*o.Hgsm + wv*o.Hvsm
			o.Hcsm = []; # free memory
			o.Hxsm = [];
		else:
			o.ham1 = g*o.Hgsm + wv*o.Hvsm
		o.Hgsm =[]; o.Hvsm = [];
	elif loopover =='wr':
		fsolve = fdiagw;# for diagonalisation
		if (detuning == 1 ):
			o.ham1 = wc*o.Hcsm + wx*o.Hxsm + wv*o.Hvsm + lambda0*wv*o.Hbsm + wv*lambda0**2*o.sft;
			o.Hcsm = []; # free memory
			o.Hxsm = [];
		else:
			o.ham1 = wv*o.Hvsm + lambda0*wv*o.Hbsm + wv*lambda0**2*o.sft;
		o.Hbsm = [];
		o.Hvsm = [];


	#lambin0 = np.linspace(lmin,lmax, nlmax);
	#print('lmin,lmax, nlmax = ',lmin,lmax, nlmax)
	#lambin0 = o.lambin0;
	lambin = []; il=-1;
	for lamb00 in lambin0:
		il=il+1
		lambin.append([il,lamb00])		
	o.lambin = lambin;
	if nlmax<Np:
		Npdiag = nlmax;
	else:
		Npdiag = Np

	
	# Diagonalisation:
	for lamb00 in lambin:
		il, evalues, o.evecs = fsolve(lamb00);
		print(' matrix elements for il = '+str(il+1)+'/'+str(nlmax)+' done ... ');

		# get shifted evalues
		if detuning==1:
			evalues = evalues - wx; # shift all E_m by E_x:
		else:
			evalues = evalues - eshft; # shift all E_m by E_x:	
		cos2th = [];

		for i in range(nstates):
			cos2th.append(getCos2th(i));


		# Herrera and Spano absorption
		#if 1: # define logical to control this and set here.
		#	pool=Pool(Np);
		#	hsmelem = pool.map(EiEjmatelem, range(nstates));
		#	pool.close();

		# combine the results for writing output file:
		absOUT = np.zeros((nstates,4));
		#print("o.Nv1l[0],o.Nv2l[0] = ",o.Nv1l[0],o.Nv2l[0])
		absOUT[:,0] = evalues;
		absOUT[:,1] = np.abs(o.evecs[0,:])**2;
		absOUT[:,2] = o.evecs[n1,:]; #hsmelem #
		absOUT[:,3] = cos2th;
	
		fabsorption = dumy+"/absorption.txt"
		f=open(fabsorption,'ab');
		if loopover =="lambda0":
			lamb0 = o.lambin0[il];
			wr0 = wr;
		else:
			wr0 = o.lambin0[il];
			lamb0 = lambda0;

		header = " "+str(nstates)+" "+str(n)+" "+str(m)+" "+str(mx);
		if (detuning==0):
			header += " "+str(wr0)+" "+str(lamb0)+" "+str(eshft)+" "+str(eshft)+" "+str(wv);
		else:
			header += " "+str(wr0)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
		np.savetxt(f, absOUT,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ', header=header,comments='#')
		f.close();
		f=open(fabsorption,'at')
		print('    ',file=f)
		print('    ',file=f)
		f.close()		



	print(" absorption data calculated for postprocessing!");
	return
#***************************************************

