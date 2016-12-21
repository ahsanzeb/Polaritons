
# waterfall plot

import numpy as np
import matplotlib.pyplot as plt

#****************************************
# for labels: paramteres for the abs spectrum we want to plot.
#****************************************
# waterfall plot
ei=-2; ef= 5; # energy window
lcol = 'g'; # line colour
fcol = 'c'; # fill colour exact
lcol2 = 'k'
fcol2 = 'y'; # fill colour analytical
trans = 0.3; # transparency factor; fill colour analytical
lzerocol = lcol2; # zero line color
lw=1
lw2 = 1;
show=1;  # 
pltGr = 0;  # plot green analytical?

fig = plt.figure(facecolor='w');
ax = fig.add_subplot(111, axisbg='w');
xline = np.linspace(ei,ef,10); yline = np.zeros(10);
#****************************************
fin = open("data/n-all/abs-vs-w-td.txt",'r');
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
	# Plot the line and fill under it: increase the z-order each time
	# so that lower lines and their fills are plotted over higher ones
	#ax.plot(a[:,0],a[:,1]+offset, lcol, lw=lw, zorder=(iy+1)*2);
	#ax.fill_between(a[:,0], a[:,1]+offset, offset, facecolor=fcol, lw=0, zorder=(iy+1)*2-1);
	#ax.plot(a[:,0],a[:,3]+offset, lcol2, lw=lw2, zorder=(iy+1)*2);
	#ax.fill_between(a[:,0], a[:,3]+offset, offset, facecolor=fcol2, lw=0, zorder=(iy+1)*2-1,alpha=trans);

	#ax.plot(xline,yline+offset, lzerocol, lw=lw2); # line to complete horizontal lines in smaller E window data
#----------------------------	
	if ldat > istop+3: #2 empty lines + 1???
		istart = istop+2; ignu += 1;
	else:
		carryon=False;
#****************************************



# max of a given quantity for all data set:
maxval = np.amax(np.array(maxvals),axis=0);

fname='data/n-all/abs-vs-t-dt-normalised.txt';
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
		if pltGr:
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


