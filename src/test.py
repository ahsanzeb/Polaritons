
# test to compare old and new modules for basis related calculations

import globalvariables as o
import basis0, basis, basisnew, mapping
import numpy as np
from memtime import memtime
import time, sys, scipy



if 0:
	n = 3; m = 5;
	for nn in range(1,n):
		ntot = int(scipy.special.binom(m+nn,nn));
		print(ntot-1)


if 0:
	x = 0; N=10000000
	for i in range(N*10):
		x  += 1;
		y = x*x;
		if i%N == 0:
			memtime(str(i));# print memory and time info

if 1:
	for n in range(20,21):
		o.m = 8
		o.mx = o.m;
		o.Np = 8
		o.n = n
		mapping.getmapp(o.n,o.m,o.Np);
		memtime('basis');# print memory and time info		
		exit()
		
		mapfulp = o.map21; o.map21 = [];
		mapping.getmap(o.n,o.m);
		r = np.sum(np.abs(np.array(mapfulp)-np.array(o.map21)))
		if r > 0:
			print(' r = ',r,'??????????')
		else:
			print('n = ',n,' r = 0; PASSED')
	exit()



if 1:
	n=14
	o.n = n;
	o.m = 8; o.mx = o.m; o.Np = 8
	o.corrtd =1; o.ld = 0.5;
	basis.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('basis');# print memory and time info		
	basisnew.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('newbasis');# print memory and time info




exit()

for n in range(14,15):
	o.n = n;
	o.m = 8; o.mx = o.m; o.Np = 8
	o.corrtd =1; o.ld = 0.5;

	memtime('rough');# print memory and time info
	basis0.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('oldbasis');# print memory and time info

	Factlist = o.Fact1l; o.Fact1l = [];
	Nvlist = o.Nv1l; o.Nv1l = [];
	Norm1,Norm2,Norm3 = o.Norm1l, o.Norm2l, o.Norm3l
	o.Norm1l, o.Norm2l, o.Norm3l = [],[],[];
	mapbasis = o.map21; o.map21 = [];
	memtime('rough');

	basisnew.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('newbasis');# print memory and time info
	mapping.getmap(o.n,o.m);
	mapp = o.map21; o.map21 = [];
	memtime('map');# print memory and time info
	mapping.getmapp(o.n,o.m,o.Np);
	memtime('mapp');# print memory and time info

	#mapping.getmap(o.n,o.m);

	r0 = np.sum(np.abs(np.array(mapbasis)-np.array(o.map21)))
	r1 = np.sum(np.abs(np.array(Factlist)-np.array(o.Fact1l)))
	r2 = np.sum(np.abs(np.array(Nvlist)-np.array(o.Nv1l)))
	r11 = np.sum(np.abs(np.array(Norm1)-np.array(o.Norm1l)))
	r12 = np.sum(np.abs(np.array(Norm2)-np.array(o.Norm2l)))
	r13 = np.sum(np.abs(np.array(Norm3)-np.array(o.Norm3l)))
	r = r0+float(r1)+r2+r12+r13+r11;
	if r>0:
		print('------    n= '+str(n)+'   ????? ... ==================????')
		print('n= '+str(n)+' r0,r1,r2,r11,r12,r13 = ',r0,r1,r2,r11,r12,r13)
	else:
		print('    ')
		print('------    n= '+str(n)+'   PASSED... ')
		print('    ')
	

