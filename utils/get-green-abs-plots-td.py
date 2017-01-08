
# if freq window is not right in the calculations,
# use this script to get FT from the correlation file
# for desired freq window/npnts.

# also generates waterfall plot

import numpy as np
import matplotlib.pyplot as plt
import decimal
import numpy as np;
from math import factorial
from numpy import pi, exp

#****************************************
# for labels: paramteres for the abs spectrum we want to plot.
#****************************************
# waterfall plot
ei=-2.5; ef= 2.5; # energy window
lcol = 'g'; # line colour
fcol = 'c'; # fill colour exact
lcol2 = 'k'
fcol2 = 'y'; # fill colour analytical
trans = 0.3; # transparency factor; fill colour analytical
lzerocol = lcol2; # zero line color
lw=1
lw2 = 1;
show=1;  # 
pltGr = 1;  # plot green analytical?


e1,e2 = ei,ef;

def getFourierTransform(Xt,ti,dd,wlist,arrow):
	Xw = []; ftnorm = 1/np.sqrt(2*pi);
	p = 1j;
	if arrow == 'B':
		p = -1j;
	for w in wlist:
		xw = 0; 	t= ti;
		for xt in Xt:
			xw += xt*exp(p*w*t);
			t += dd;
		Xw.append(xw*dd*ftnorm);
	Xw = np.array(Xw);
	return Xw
def getFrankCondonEtc(l0):
	res=[]; 	efac=[];
	if abs(l0) < 1e-8:
		for i in range(mx+1):
			if i==0: x=1;
			else: x=0;
			res.append(x);
			efac.append(i*wv -1j*gamma/2);
		return res,efac
	cont=decimal.Context(prec=15, Emax=999, clamp=1);
	decl0 = decimal.Decimal.from_float(l0);
	for i in range(mx+1):
		#x = np.exp(-l0**2)*l0**(2*i)/factorial(i);
		n0 = decimal.Decimal.from_float(factorial(i));
		fc0 = cont.power(decl0,2*i);
		x = np.exp(-l0**2)*float(fc0/n0);
		res.append(x);
		efac.append(i*wv -wv*l0**2 -1j*gamma/2); # -wv*l0**2 polaron transform term not included in JK notes
	return res,efac
# ---------------------
def Green(wlist,l0,wr):
	# ------------------
	# getFrankCondonEtc()
	FC, efac = getFrankCondonEtc(l0);# Frank-Condon factors etc
	# ------------------
	wr2 = wr*wr; delta = wc-wx;
	GRlist = [];
	for w in wlist:
		SelfE = 0;
		for i in range(mx+1):
			SelfE -= wr2*FC[i]/(w - efac[i]);
		GR = 1/(w +1j*kappa/2 - delta + SelfE);
		GRlist.append(GR);
	GRlist = np.array(GRlist);
	return GRlist
# -------------------------------------------
def fwriteFT(il,wlist, Gw, GR, fnametd, ntmax):
	# save FT of correlation
	absOUT = np.zeros((nwmax,4));
	absOUT[:,0] = wlist;
	absOUT[:,1] = np.real(Gw);
	# 2*kl*np.imag(absws)**2 + kappa*np.abs(absws)**2;
	absOUT[:,2] = 2*kl*np.real(Gw)**2 + kappa*np.abs(Gw)**2;
	absOUT[:,3] = -np.imag(GR); # Green from analytical
	f=open(fnametd,'ab');
	lamb0 = lamblist[il];
	wr0 = wrlist[il];
	header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nwmax)+" "+str(ntmax);
	header += " "+str(wr0)+" "+str(lamb0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	np.savetxt(f, absOUT,fmt='%15.10f %15.10f %15.10f %15.10f', delimiter=' ', header=header,comments='#')
	f.close();
	f=open(fnametd,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close()		
	return None

def fcorrft(il,corr):
	E1, E2 = e1+wx, e2+wx;	
	#print('using: e1,e2 = ',e1,e2)
	wlist = np.linspace(E1, E2,nwmax);
	Gw = getFourierTransform(corr,0,dt,wlist,'F');
	# green function based results:
	# wlist = wlist - wx ; # shift freq axis
	lamb0 = lamblist[il];
	wlist = np.linspace(e1,e2,nwmax);
	wr = wrlist[il];
	GR = Green(wlist,lamb0,wr); 	# analytical Green function
	return Gw, GR, wlist


# ------------------------------
fin = open("corr-vs-t-td.txt",'r');
lines = fin.readlines();
ldat = len(lines);
carryon=True;
istart = 0; ignu=0; results = []; headers = []; maxvals = [];
lamblist=[]; wrlist=[];

while carryon==True:
	# read parameters:
	line = lines[istart];
	linesplit = line.split();
	z = [int(x) for x in linesplit[1:6]];
	# n,m,mx,nw,ntmax = z[0],z[1],z[2],z[3],z[4];
	n,m,mx,nw,ntmax = z;
	z = [float(x) for x in linesplit[6:15]];
	# wr,l0,wc,wx,wv,gamma,kappa,tf,dt = z[0],z[1],z[2],z[3],z[4],z[5],z[6],z[7],z[8]
	wr,l0,wc,wx,wv,gamma,kappa,tf,dt = z;
	lamblist.append(l0);
	wrlist.append(wr); kl=kappa/2; nwmax=nw;
	# header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nw)+" "+str(ntmax);
	# header += " "+str(wr)+" "+str(l0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	# header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	# headers.append(header);
	headers.append(line);

	# kr =kappa/2; kl=kappa/2;
	wr2 = wr**2; delta = wc-wx;
	istop = istart + ntmax +1;
	a=[];
	for line in lines[istart+1:istop]:
		a.append( [ float (x) for x in line.split() ] );
	ntot = len(a);
	a=np.array(a);
	results.append(a);
	#.............................
	iy = ignu;
	offset = -iy # (ny-iy)*5
	if ldat > istop+3: #2 empty lines + 1???
		istart = istop+2; ignu += 1;
	else:
		carryon=False;
# ------------------------------

print(' read--------')

# do everything!
il = 0;
for x in results:
	corr = x[:,1] + 1j*x[:,2];
	Gw, GR, wlist = fcorrft(il,corr);
	# 1/2.5 makes GR peaks equal to the exact for lam=0 case.
	fwriteFT(il, wlist, Gw, GR/2.5, 'abs-vs-w-td-1.txt', ntmax); # write FT file
	il += 1;
# ------------------------------


fig = plt.figure(facecolor='w');
ax = fig.add_subplot(111, axisbg='w');
xline = np.linspace(ei,ef,10); yline = np.zeros(10);
#****************************************
fin = open("abs-vs-w-td-1.txt",'r');
lines = fin.readlines();
ldat = len(lines);
carryon=True;
istart = 0; ignu=0; results = []; headers = []; maxvals = [];
while carryon==True:
	# read parameters:
	line = lines[istart];
	linesplit = line.split();
	z = [int(x) for x in linesplit[1:6]];
	# n,m,mx,nw,ntmax = z[0],z[1],z[2],z[3],z[4];
	n,m,mx,nw,ntmax = z;
	z = [float(x) for x in linesplit[6:15]];
	# wr,l0,wc,wx,wv,gamma,kappa,tf,dt = z[0],z[1],z[2],z[3],z[4],z[5],z[6],z[7],z[8]
	wr,l0,wc,wx,wv,gamma,kappa,tf,dt = z;
	# header =" "+str(n)+" "+str(m)+" "+str(mx)+" "+str(nw)+" "+str(ntmax);
	# header += " "+str(wr)+" "+str(l0)+" "+str(wc)+" "+str(wx)+" "+str(wv);
	# header += " "+str(gamma)+ " "+str(kappa)+ " "+str(tf)+ " "+str(dt);
	# headers.append(header);
	headers.append(line);
	
	# kr =kappa/2; kl=kappa/2;
	wr2 = wr**2; delta = wc-wx;
	istop = istart + nw +1;
	a=[];
	for line in lines[istart+1:istop]:
		a.append( [ float (x) for x in line.split() ] );
	ntot = len(a);
	a=np.array(a);
	results.append(a);
		# normalise all data by max peak value
	maxval = np.amax(a[:,1:],axis=0);
	# print("maxval = ",maxval);
	maxvals.append(maxval);
#.............................
	iy = ignu;
	offset = -iy # (ny-iy)*5
#----------------------------	
	if ldat > istop+3: #2 empty lines + 1???
		istart = istop+2; ignu += 1;
	else:
		carryon=False;
#****************************************



# max of a given quantity for all data set:
maxval = np.amax(np.array(maxvals),axis=0);

fname='abs-vs-t-dt-normalised.txt';
i = 0;
for a in results:
	a[:,1] = a[:,1]/maxval[0];
	a[:,2] = a[:,2]/maxval[1];
	a[:,3] = a[:,3]/maxval[2];
	header = headers[i];
	f=open(fname,'ab');
	np.savetxt(f, a,fmt='%15.10f%15.10f%15.10f%15.10f', delimiter=' ',header=header,comments=' ')
	f.close();
	f=open(fname,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close();
	# waterfall plot ,comments='#'
	iy = i
	if show:
		offset = -iy # (ny-iy)*5
		w = a[:,0];
		ax.plot(w,a[:,1]+offset, lcol, lw=1, zorder=(iy+1)*2);
		ax.fill_between(w, a[:,1]+offset, offset, facecolor=fcol, lw=0, zorder=(iy+1)*2-1);
		if pltGr and i==len(results)-1:
			ax.plot(a[:,0],a[:,3]+offset, lcol2, lw=lw2, zorder=(iy+1)*2);
			ax.plot(xline,yline+offset, lcol, lw=1);
	i += 1;

if show:
	ttl = "${HTC\,Model:\,\,N= "+str(n)+"\,\kappa= "+str(kappa)+"\,\gamma= "+str(gamma)+"}$";
	plt.xlim(ei,ef);
	plt.ylim(-ignu-0.2,1.2);
	plt.xlabel("$\omega-\omega_0$",color='k',fontsize=24);
	plt.ylabel("$\Im G^R$",color='k',fontsize=24);
	plt.title(ttl,color='k',fontsize=16);
	plt.yticks([], [])
	plt.tight_layout();
	plt.savefig('abs-from-td-waterfall.pdf', format='pdf');
	plt.show();


