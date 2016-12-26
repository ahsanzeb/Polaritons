
import globalvariables as o
import basisnew
import basis0
import numpy as np
from memtime import memtime
import time, sys





o.m, o.mx = 5,5;
o.Np = 8;

o.corrtd =0; o.ld = 0.5;

for n in range(10,12):
	o.n = n;

	basis0.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('basis old');# print memory and time info

	Factlist = o.Fact1l;
	Nvlist = o.Nv1l
	Norm1,Norm2,Norm3 = o.Norm1l, o.Norm2l, o.Norm3l

	# memtime('----');# print memory and time info

	basisnew.fbasis(o.n,o.m,o.mx,o.Np);
	memtime('basis new');# print memory and time info

	r1 = np.sum(np.abs(np.array(Factlist)-np.array(o.Fact1l)))
	r2 = np.sum(np.abs(np.array(Nvlist)-np.array(o.Nv1l)))
	r11 = np.sum(np.abs(np.array(Norm1)-np.array(o.Norm1l)))
	r12 = np.sum(np.abs(np.array(Norm2)-np.array(o.Norm2l)))
	r13 = np.sum(np.abs(np.array(Norm3)-np.array(o.Norm3l)))
	r = r1+r2+r12+r13+r11;
	if r>0:
		print('------    n= '+str(n)+'   ????? ... ==================????')
		print('n= '+str(n)+' r1,r2,r11,r12,r13 = ',r1,r2,r11,r12,r13)
	else:
		print('    ')
		print('------    n= '+str(n)+'   PASSED... ')
		print('    ')
	

