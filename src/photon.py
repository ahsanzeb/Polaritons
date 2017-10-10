
import globalvariables as o
import numpy as np

# define local to avoid writing o. every time!
n,m,mx,Np = o.n, o.m, o.mx,o.Np
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv
dumy = o.dumy;
n1,n2 = o.n1,o.n2;
n1fsym = o.n1fsym;
nlmax = o.nlmax;
loopover= o.loopover;
lambin0 = o.lambin0;
lambda0 = o.lambda0;

listn2 = o.listn2;
listn3 = o.listn3;
# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------




def photfrac():
	print(' ====> calculating photon fractions... ')
	# -----------------------------
	if loopover == 'lambda0':
	# parameters:
		# param = np.array([n,m,wv,wr,lamb]);
		param = np.array([n,m,wv,wr])
	# descriptive output file names
		fout0=dumy+"/fractions-N="+str(n)+"-M="+str(m)+"-wr="+str('%.1f'%(wr))+".dat"
	else: # if loopover == 'wr':
	# parameters:
		# param = np.array([n,m,wv,wr,lamb]);
		param = np.array([n,m,wv,lambda0])
	# descriptive output file names
		fout0=dumy+"/fractions-N="+str(n)+"-M="+str(m)+"-l0="+str('%.1f'%(lambda0))+".dat"
	# -----------------------------

	dm0 = np.zeros(nlmax);
	for il, evalu, evec in o.eigvv:
		dm0[il] = np.sum(evec[0:n1][0]);	
	writephot(dm0,param,lambin0,fout0);
	return

#-------------------------------------
# for writing output files 
#-------------------------------------
def writephot(dm,param,lambin0,fout):
	f=open(fout,'ab')
	np.savetxt(f,param[None],  fmt='%5i, %5i, %10.6f, %10.6f')
	for il in range(nlmax):
		x = dm[il];
		l0 = lambin0[il];
		res = np.array([l0,x,(1-x)*1.0/n,(1-x)*(n-1)*1.0/n])
		np.savetxt(f,res[None], fmt='%15.10f', delimiter='')
	f.close()
	return
#-------------------------------------

