

import globalvariables as o
import numpy as np

# define local to avoid writing o. every time!
dumy = o.dumy;
loopover= o.loopover;
dumy = o.dumy;
lambin = o.lambin;
justenergy = o.justenergy;
nlmax = o.nlmax;

# add as comments in output file
n,m,mx = o.n, o.m, o.mx
wr,wx,wc,wv = o.wr, o.wx, o.wc, o.wv

# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------

print(' ====> writing energy output file ');

if loopover == "lambda0":
	fenergy = dumy+"/energy-vs-lambda0.txt"
else:
	fenergy = dumy+"/energy-vs-wr.txt"
	

# energy = [[0.0]*2 for il in range(nlmax)];
energy=[];
if justenergy:
	for il, evalu in o.eigvv:
		# energy[il][0] = lambin[il][1]
		# energy[il][1] = evalu
		energy.append([lambin[il][1]]+evalu.tolist())
else:
	energy = [[0.0]*2 for il in range(nlmax)];
	for il, evalu, evec in o.eigvv:
		energy[il][0] = lambin[il][1]
		energy[il][1] = evalu;
energy = np.array(energy)
# writing energy file 
f=open(fenergy,'ab')
np.savetxt(f, energy)
# np.savetxt(f, energy, fmt='%10.7f %15.10f')
f.close()
f=open(fenergy,'at')
print('    ',file=f)
print('    ',file=f)
print('    ',file=f)
# print('# # # Comment lines for gnuplot index',file=f)
# print('# # # '+" N="+str(n)+", M="+str(m),file=f)
# print('# # # # # # # # # # # # ',file=f)
f.close()


