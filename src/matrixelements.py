


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
detuning, lamb0, eshft, nstates = o.detuning, o.lamb0, o.eshft, o.nstates

lmin,lmax, nlmax = o.lmin, o.lmax,  o.nlmax;
loopover= o.loopover;
lambda0= o.lambda0;


# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

#-----------------------------------------------
#ROUGH COMMENTS
#-----------------------------------------------
	# "evec0hc: only exciton block elements needed"
	# "so translate the index back by n1 when using j"

	# "in main prog: mk pool over eigvv (eigvv is a list), apply to absorption()"
	# "arg: also give evec0hc"
	# "FINAL: given indeces to map, use the same pool that knows eigvv"
	#for il, evalu, evec in eigvv:
		#evec # used in: phot_fraction(), abs_matelem()
#-----------------------------------------------


#-----------------------------------------------
# These variable should be define before calling MatElemAndPhotFrac():
# evec0hc, eigvv, (and obviously, n1, n2, m, o.map21)
#-----------------------------------------------
# use in pool map
#def	MatElemAndPhotFrac(ierange):
# Input: ierange: [ie1,ie2];
#	listout = []; pf = [];
#	ie1 = ierange[0]; ie2 = ierange[1];
#	for i in range(ie1,ie2):
#		evec = o.evecs[i];
#		ma0 = abs_matelem(i);
#		pf = phot_fraction(i); 
#		listout.append([ma0,pf]); 
#	return listout;
#--------------------------


def getCos2th(i):
	cos2th = 0;
	# first n1 basis have a photon
	for p in range(n1): 
		cos2th += abs(o.evecs[p,i])**2
	return cos2th

# matrix elements between LP and other eigenstates of 1 excitation space (the one we calculate in the code)
def abs_matelem(i):# uses eigestate0hc[n1:]
	matelem = 0.0;
	evec = o.evecs[i];
	for jj in range(n2):
		for mj in range(0,m+1):
			# n1+ mj*n2 + jj translated by -n1; exciton block basis index
			j= mj*n2 + jj; 
			ii = o.map21[jj,mj]; # photon block basis index
			x = evec[ii]*evec0hc[j];
			matelem += x  #
			# activate for franck condon factors comparison
			# if o.Nv2l[jj]==0:
			# 	if x > 1e-15:
			# 		print("mj, o.Nv1l[ii], evec[ii] = ",mj, o.Nv1l[ii], evec[ii])
		# extra basis: diag in mj,jj
		for mj in range(m+1,mx+1): 
			# coor in exciton block
			j = mj*n2 + jj; 
			# coor in photon block
			ii = n1fsym +(mj-m-1)*n2 + jj;
			x = evec[ii]*evec0hc[j];
			matelem +=  x
	return abs(matelem)**2
#--------------------------
def phot_fraction(i):
	cos2th = 0;
	#evec = o.evecs[i];
	# first n1 basis have a photon
	for p in range(n1): 
		cos2th += abs(o.evecs[i][p])**2
	return cos2th
#--------------------------
















# -----------------------------------------------------
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


	lambin0 = np.linspace(lmin,lmax, nlmax);
	print('lmin,lmax, nlmax = ',lmin,lmax, nlmax)
	o.lambin0 = lambin0;
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
			#print("wx = ",wx)
			evalues = evalues - wx; # shift all E_m by E_x:
		else:
			#print("wx = ",eshft)
			evalues = evalues - eshft; # shift all E_m by E_x:	
		pool=Pool(Np); # refresh pool
		# photon fraction:
		cos2th = pool.map(getCos2th, range(nstates))
		# coeff of eigenstates in photon and exciton sectors with 0 phonon: basis indeces = 0,n1
		# Aphot = o.evecs[0,:]; # all eigenstates, columns
		# Aexc = o.evecs[n1,:];
		# combine the results for writing output file:
		absOUT = np.zeros((nstates,4));
		#print("o.Nv1l[0],o.Nv2l[0] = ",o.Nv1l[0],o.Nv2l[0])
		absOUT[:,0] = evalues;
		absOUT[:,1] = o.evecs[0,:];
		absOUT[:,2] = o.evecs[n1,:];
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
		pool.close(); # close this pool

	print(" absorption data calculated for postprocessing!");
	return
#***************************************************

