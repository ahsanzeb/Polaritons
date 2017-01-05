
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import math

# ---------------------# ---------------------
# Input:
# ---------------------# ---------------------
# decay rates for cavity and excitons:
kr= 0.05;
gamma = 0.05;
kl = kr; 

# bare excitonic dipole transition matrix elemet
alpha = 0.0; # depends on probe beam intensity, I think!
# I think, should be less than a few % of rabi frequency.
# should not it be much bigger than gamma, because that gamma is due to coupling of atom with just free space radiation modes but this alpha is due to coupling with resonant probe beam?!

kappa = 0.5*(kr+kl);


dw = 0.001; # frequency bin width:
wi = -2.5;
wf = 2.5;# set to -ve to include all eigenvectors

# waterfall plot
ei=-1.5; ef= 1.5; # energy window
lcol = 'b'; fcol = 'b'; # line and filling colours

fig ='false';
writeout = 'true';
Np= 8; # no. of processes for pool
#import param
#n=param.n # n = no. of sites
#m=param.m # m = max no. of phonons per site
#l0=param.lmin # lambda_0 

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbg='w')
xline = np.linspace(ei,ef,10);
yline = np.zeros(10);


fabsorption = "absorption.txt";
fname = "abs-vs-w.txt";
#fabsorption = "data/n-"+str(n)+"/absorption.txt"
#fname = "data/n-"+str(n)+"/abs-vs-w.txt";
# ---------------------# ---------------------
# Functions
# ---------------------# ---------------------
# flattens output of pool.map to make a list of elements
flat = lambda l: [item for sublist in l for item in sublist];
# ---------------------
# range for floats:
def frange(x, y, dx):
	out=[];
	while x < y:
		out.append(x);
		x += dx;
	return np.array(out)
# ---------------------
def denom(i):
	cos2th = a[i,3]; # photon fraction
	gm = (kr+kl)*cos2th + gamma*(1-cos2th);
	return a[i,0] + 1j*gm*0.5;
# ---------------------
def absw(w):
	res = 0.0;
	for i in range(ntot):
		res += abs((a[i,1] + alpha*a[i,2]))**2/(w - EGamma[i]);
	return res
# ---------------------
def getFrankCondonEtc():
	res=[]; 	efac=[];
	for i in range(5*mx+1):
		x = np.exp(-l0**2)*l0**(2*i)/math.factorial(i);
		res.append(x)
		efac.append(i*wv -wv*l0**2 -1j*gamma/2); # -wv*l0**2 polaron transform term not included in JK notes
	return res,efac
# ---------------------
def Green(w):
	SelfE = 0;
	for i in range(5*mx+1):
		SelfE -= wr2*FC[i]/(w - efac[i]);
	GR = 1/(w +1j*kappa/2 - delta + SelfE);
	return GR
# ---------------------






# ---------------------
# calculations:
# ---------------------
# read "absorption.txt"
# Format: Em-E0, |<0|a|m>|^2, cos2th_m
# fabsorption = "absorption.txt";
	              
fgnu=open("gnuplot.txt",'wt'); gline='plot ';

fin = open(fabsorption,'r');
lines = fin.readlines();
ldat = len(lines);
carryon=True;
istart = 0; ignu=0;
results = []; maxvals = []; headers = [];
while carryon==True:
	print(ignu)
	# read parameters:
	line = lines[istart];
	linesplit = line.split();
	z = [int(x) for x in linesplit[1:5]];
	# print(len(z))
	# print(z);
	nstates,n,m,mx = z[0],z[1],z[2],z[3];
	z = [float(x) for x in linesplit[5:10]];
	# print(len(z))
	# print(z);
	wr,l0,wc,wx,wv = z[0],z[1],z[2],z[3],z[4];
	wr2 = wr**2; delta = wc-wx;
	header = line;
	# header = " "+str(nstates)+" "+str(n)+" "+str(m)+" "+str(mx);
	# header += " "+str(wr)+" "+str(l0)+" "+str(wc)+" "+str(wx)+" "+str(wv)+" "+str(nw);

	istop = istart + nstates +1;
	# print("strt,istop,nstates=",istart+1, istop,nstates)
	a=[];
	for line in lines[istart+1:istop]:
		a.append( [ float (x) for x in line.split() ] );
	ntot = len(a);
	a=np.array(a);
	#----------------------------
	#calcabs();
	#----------------------------
	if wi>=0 or wi<=a[0,0]:
		w1 = a[0,0] - 10*(kappa+gamma); # min value of Em-Ex
	else:
		w1 = wi;
	if wf<=0 or wf>=a[ntot-1,0]:
		w2 = a[ntot-1,0] + 10*(kappa+gamma);# max value of Em-Ex
	else:
		w2 = wf;		
	# frequency values for which absorption is to be calcualted
	w = frange(w1, w2, dw);
	nw = len(w);
	headers.append(header+" "+str(nw));
	# calculate absorption spectrum:
	# if not efficient, call functions less by chunking the argument list
	pool=Pool(Np) # Np processes
	# denominator terms: list of (E_m-E_0)+I*gamma_m
	EGamma = pool.map(denom,range(ntot));
	#EGamma = flat(EGamma);
	# green function based results:
	FC, efac = getFrankCondonEtc(); 	# Frank-Condon factors etc
	# refresh pool
	pool.close();
	pool=Pool(Np)
	# "absorption" for all frequencies:
	absws = pool.map(absw,w); # Green from Numerics
	# -------------------------------
	GR = pool.map(Green,w);
	# -------------------------------
	# ---- write output file -----
	# E, ImG, A(w), ImGa
	result = np.zeros((nw,4));
	result[:,0] = w;
	result[:,1] = np.imag(absws);
	result[:,2] = 2*kl*np.imag(absws)**2 + kappa*np.abs(absws)**2;
	result[:,3] = -np.imag(GR); # Green from analytical
	# normalise all data by max peak value
	maxval = np.amax(result[:,1:],axis=0);
	# print("maxval = ",maxval);
	maxvals.append(maxval);
	results.append(result);

	# flam=open("lam-"+str(l0)+".txt",'ab')
	# np.savetxt(flam, result,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ')
	# flam.close();
	pool.close();
	slabel = "G(w): n,m,mx="+str(n)+", "+str(m)+", "+str(mx)
	gline += ' "abs-vs-w.txt" i '+str(ignu);
	gline += ' w l lw 3 t '+repr(slabel)
	gline += ', "abs-vs-w.txt" i '+str(ignu);
	gline += ' u 1:3 w l lw 2 t "A(w)"'
	gline += ', "abs-vs-w.txt" i '+str(ignu);
	gline += ' u 1:4 w l lw 1 t "MEx:G(w)"'
	#----------------------------	
	if ldat > istop+3: #2 empty lines + 1???
		istart = istop+2; ignu += 1;
		gline +=', '
	else:
		carryon=False;
	#print(istop, ldat)
#print("********")
print(gline,file=fgnu)
print('    ',file=fgnu)
fgnu.close()		




# max of a given quantity for all data set:
maxval = np.amax(np.array(maxvals),axis=0);
print("maxval = ",maxval);

i = 0;
for result in results:
	result[:,1] = result[:,1]/maxval[0];
	result[:,2] = result[:,2]/maxval[1];
	result[:,3] = result[:,3]/maxval[2];
	header = headers[i];
	f=open(fname,'ab');
	np.savetxt(f, result,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ',header=header,comments=' ')
	f.close();
	f=open(fname,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close();
	# waterfall plot
	iy = i
	offset = -iy # (ny-iy)*5
	w = result[:,0];
	ax.plot(w,result[:,1]+offset, lcol, lw=1, zorder=(iy+1)*2);
	ax.fill_between(w, result[:,1]+offset, offset, facecolor=fcol, lw=0, zorder=(iy+1)*2-1);
	ax.plot(xline,yline+offset, lcol, lw=1);
	i += 1;

print("Absorption spectrum calculated! .... ")

ttl = "${HTC Model: Absorption N= "+str(n)+" M= "+str(m)+" \kappa_{L,R}= "+str(kr)+" \gamma= "+str(gamma)+"}$";

plt.xlim(ei,ef);
plt.ylim(-ignu-0.2,1.2);
plt.xlabel("$\omega-\omega_0$",color='k',fontsize=24);
plt.ylabel("$\Im G^R$",color='k',fontsize=24);
plt.title(ttl,color='k',fontsize=16);
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
plt.tight_layout();
plt.savefig('fig-waterfall.pdf', format='pdf');
plt.show()

