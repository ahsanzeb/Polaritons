

from scipy.special import binom
import numpy as np

maxsize = 10e6;
m,mx = 3,3;
for n in np.linspace(5,100,20):
	n1 = binom(m+n, m);
	n2 = binom(m+n-1, m);
	ntot = n1 + (mx+1)*n2;	
	x=int(n), int(ntot)
	print(x);

exit()


n=30
if 1:
	for m in range(1,mx+1):
		n1fsym = binom(m+n, m);
		n2 = binom(m+n-1, m);
		n1 = n1fsym + n2*(mx-m);
		ntot = n1 + (mx+1)*n2;	
		x=int(m), int(n1), int(n2*(mx+1)), int(ntot)
		print(x);

exit()



for n in np.linspace(5,100,21):
	print(' ------ n = ',n,' -------')
	for m in range(1,mx+1):
		n1fsym = binom(m+n, m);
		n2 = binom(m+n-1, m);
		n1 = n1fsym + n2*(mx-m);
		ntot = n1 + (mx+1)*n2;	
		if ntot <= maxsize:
			mm = m;
			x=int(m), int(n1), int(n2*(mx+1)), int(ntot)
	print(x);
	print('Full Hilbert: ',(m+1)**n)

