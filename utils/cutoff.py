

from scipy.special import binom
import numpy as np
from math import factorial as fct
import decimal
import matplotlib.pyplot as plt
from numpy import exp 


def getmlist(mcut):
	mlist = [];
	cont=decimal.Context(prec=15, Emax=999, clamp=1);
	for lam0 in lamlist:
		decl0 = decimal.Decimal.from_float(lam0);
		occups=[];
		for i in range(mmax):
			xa = cont.power(decl0,i);
			y = decimal.Decimal(fct(i));
			z = xa/y;
			x = float(z) * exp(-lam0**2);
			if i < lam0**2 + lam0:
				occups.append([i,x]);
			elif x > mcut:
				occups.append([i,x]);
		mlist.append(occups[-1][0]);
	return np.array(mlist);


lamlist = np.linspace(0.1,2.5,10)
mmax = 30;
cont=decimal.Context(prec=15, Emax=999, clamp=1);

fig = plt.figure();
ax= fig.gca();
ax.set_yscale('log',basey=10)
for lam0 in lamlist:
	decl0 = decimal.Decimal.from_float(lam0);
	occups=[];
	for i in range(mmax):
		xa = cont.power(decl0,i);
		y = decimal.Decimal(fct(i));
		z = xa/y;
		x = float(z) * exp(-lam0**2);
		if i < lam0**2 + lam0:
			occups.append([i,x]);
		elif x > 10e-15:
			occups.append([i,x]);
	occups = np.array(occups);
	ax.plot(occups[:,0],occups[:,1],label=str(np.round(lam0,3)))

#plt.legend();
#plt.show()

fig = plt.figure();
mldat = []; dat = [];
for j in range(5,16,2):
	valcut = 10**(-j);
	mlist = getmlist(valcut);
	plt.plot(lamlist, mlist,'-o',label="${10^{-"+str(j)+"}}$")
plt.title("basis cutoff for different tail values of  Poisson disptribution")
plt.legend(loc=2);
plt.show()





