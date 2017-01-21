
import scipy
import os
path = os.getcwd();

def 	prntmsg(n,m,mx,niter,lnlist):
	print('');
	print('======= n = '+str(n)+' ==== '+str(niter+1)+'/'+str(lnlist)+' ======== ');
	print('        (m, mx = '+str(m)+', '+str(mx)+')')
	print('');
	return

def createoutdir(N,diffoutdir):
	# create different dir for output?
	if diffoutdir:
		dumy="data/n-"+str(N)
		dpath = path+"/"+dumy
		os.makedirs(dpath, exist_ok=True)
	else:
		dumy="data/n-all"
		dpath = path+"/"+dumy
		os.makedirs(dpath, exist_ok=True)
	return dumy

def setNp(n,m,Np):
	n3 = int(scipy.special.binom(m+n-2,m));
	if n3 < Np:
		print(' n3 < Np ====> o.Np=n3');
		Np = n3;
	return Np


