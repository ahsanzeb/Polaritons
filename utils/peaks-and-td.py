
# waterfall plot

import numpy as np
import matplotlib.pyplot as plt
#****************************************
# for labels: paramteres for the abs spectrum we want to plot.
#****************************************
# waterfall plot
ei=-1; ef= 1; # energy window
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

getpeaks = 1; # get peak positons and show in plot.


def findpeaks(Y):
	cand = []; increasing = 1;
	y1 = 0; i = 0;
	for y in Y[1:]:
		i += 1;
		if y > y1: increasing = 1;
		if increasing and y < y1:
			cand.append(i-1)	; # i or i-1?
			increasing = 0;
		y1 = y;
	return cand	


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
	if ignu>0: results.append(a);
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

peakdata = [];

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

		if getpeaks:
			Y = a[:,1];
			#peaks = getpeakpositions(Y);
			peaks = findpeaks(Y);
			peakdata.append([w[peaks],Y[peaks]])
			#plt.plot(w,Y+offset)
			plt.scatter(w[peaks], Y[peaks]+offset, marker='x', color='r', s=40)
			#ax.plot(w,Y, lcol, lw=1, zorder=(iy+1)*2);
			#plt.xlim(-1.2,1.2);
			#plt.show();
			#exit()

		ax.plot(w,a[:,1]+offset, lcol, lw=1, zorder=(iy+1)*2);
		ax.fill_between(w, a[:,1]+offset, offset, facecolor=fcol, lw=0, zorder=(iy+1)*2-1);
		if pltGr and i==len(results)-1:
			ax.plot(a[:,0],a[:,3]+offset, lcol2, lw=lw2, zorder=(iy+1)*2);
			ax.plot(xline,yline+offset, lcol, lw=1);
			Y = a[:,3];
			peaks = findpeaks(Y);
			peakdata.append([w[peaks],Y[peaks]])
			#plt.plot(w,Y+offset)
			plt.scatter(w[peaks], Y[peaks]+offset, marker='x', color='r', s=40)

			
	i += 1;




if show:
	ttl = "${HTC\,Model:\,\kappa= "+str(kappa)+"\,\gamma= "+str(gamma)+"}$";
	plt.xlim(ei,ef);
	plt.ylim(-ignu-0.2+1,1.2);
	plt.xlabel("$\omega-\omega_0$",color='k',fontsize=24);
	plt.ylabel("$\Im G^R$",color='k',fontsize=24);
	plt.title(ttl,color='k',fontsize=16);
	plt.yticks([], [])
	plt.tight_layout();
	plt.savefig('abs-from-td-waterfall.pdf', format='pdf');
	plt.show();


if getpeaks:

	# find max size of w (no. of peaks)
	lm = 0; i= 0; lmin=len(peakdata[0][0]);
	for w,Y in peakdata:
		lm = max(lm,len(w));
		lmin = min(lmin,len(w))
		i += 1;
	im = i;
	pos = np.zeros((im+1,lm));
	val = np.zeros((im+1,lm));
	track = np.zeros((im+1,lmin));
	# make np array
	i= 0
	for w,Y in peakdata:
		l = len(w);
		pos[i,0:l] = w
		val[i,0:l] = Y
		i += 1;

	nlist=[2,3,4,5,6,8,10,12,14,16,18,20,22];
	iGr = len(nlist)-1;
	alpha = 0.1; alphaG = 0.1;
	datindlist = []; peaklst = [];
	w0 = peakdata[0][0];
	for pk in w0:
		i = 1; datind = []; peakpos = [];
		pk1 = pk;
		for w,Y in peakdata[1:]:
			for pk2 in w:
				pks = [pk>ei,pk1>ei,pk2>ei,pk<ef,pk1<ef,pk2<ef];
				shft = alpha;
				if i==iGr:
					shft = alphaG;
				if all(pks) and abs(pk1-pk2) < shft:
				#if pk > ei and pk < ef and pk1 > ei and pk2 < ef and pk2 > ei and pk1 < ef and abs(pk1-pk2) < alpha:
					#print(pk, i, pk1, pk2)
					datind.append(nlist[i]);
					peakpos.append(pk2);
					pk1 = pk2;
			i += 1;
		#print(datind)
		if len(datind) > 0:
			datindlist.append([datind,peakpos])
			#peaklst.append(peakpos)


	fig = plt.figure()
	ax = plt.gca();
	ax.set_xscale('log')
	ax.set_yscale('log')

	ygrd = [];
	plt.ylim(3.1,5);
	for datind,peakpos in datindlist:
		ax.plot(1/np.array(datind), np.array(peakpos)+4, 'or',ms=8);
		ax.plot(1/np.array(datind[-1]), np.array(peakpos[-1])+4, 'sg',ms=8);
		ygrd.append(np.array(peakpos[-2])+4)
	ax.set_yticks(ygrd); #np.arange(0,10,0.1));
	ax.set_yticklabels(np.round(np.array(ygrd)-4,2));	
	xticss = []; xlabels=[];
	for n in nlist[1:-1][::-1]:
		xlabels.append(str(n));
		xticss.append(1/n);
	xlabels.append('G');
	xticss.append(1/nlist[-1]);

	ax.set_xticks(xticss)
	ax.set_xticklabels(xlabels)
	plt.xlim(1/(nlist[-1]+1),1/5.5);
	ax.yaxis.grid(color='r');


	ttl = "Log-Log plot of peak position vs 1/N";
	plt.xlabel("$N$",color='k',fontsize=15);
	plt.ylabel("Peaks' positions",color='k',fontsize=15);
	plt.title(ttl,color='k',fontsize=16);
	plt.tight_layout();
	plt.savefig('abs-peaks-convergence.pdf', format='pdf');
	plt.show();

	exit()



