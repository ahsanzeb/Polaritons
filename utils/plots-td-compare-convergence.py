
# waterfall plot

import numpy as np
import matplotlib.pyplot as plt

#****************************************
# for labels: paramteres for the abs spectrum we want to plot.
#****************************************
# waterfall plot
ei=-1.5; ef= 2.5; # energy window
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
nm = len(results)//2;
b00 = results[0][:,1];
b05 = results[nm][:,1];
nw = len(b00);
resld00=[];
resld05=[];
a00 = results[1][:,1];
a05 = results[nm+1][:,1];
res00=[];res05=[];
for a in results:
	if 1:
		a[:,1] = a[:,1]/maxval[0];
		a[:,2] = a[:,2]/maxval[1];
		a[:,3] = a[:,3]/maxval[2];


	if i>0 and i <nm:
		resld00.append(np.sqrt(np.sum((a[:,1] - b00)**2))/nw);
		res00.append(np.sqrt(np.sum((a[:,1] - a00)**2))/nw);
		a00 = a[:,1];
	elif i >nm:
		resld05.append(np.sqrt(np.sum((a[:,1] - b05)**2))/nw);	
		res05.append(np.sqrt(np.sum((a[:,1] - a05)**2))/nw);
		a05 = a[:,1];
	i +=1


plt.xlim(2,30)
plt.plot(range(1,nm),res00,'sg',ms=8,label='ld = 0.0');
plt.plot(range(1,nm),res05,'or',ms=8,label='ld = 0.5');
plt.show()


if 0:
	plt.figure()
	plt.semilogy(range(1,nm),res00,'sg',ms=8,label='ld = 0.0');
	plt.semilogy(range(1,nm),res05,'or',ms=8,label='ld = 0.5');

plt.semilogy(range(1,nm),resld00,'sy',ms=8,label='ld = 0.0');
plt.semilogy(range(1,nm),resld05,'ob',ms=8,label='ld = 0.5');

plt.show()
exit()

plt.figure()

plt.xlim(1,50);
plt.xlabel('Vib cutoff $M$',fontsize=20);
txt = '$Residual(M):=\sqrt{\sum_{\omega_i}{|G(\omega_i)_M-G(\omega_i)_{50}|^2}}/\sum_{\omega_i} 1$';
txt2 = r'$Residual(M):=\sqrt{\overline{ (G_M(\omega)-G_{50}(\omega) )^2 }}$';
plt.title('Convergenvce test for absorption \n'+txt,fontsize=18);

plt.ylabel('$Residual(M)$',fontsize=20);
plt.legend();
plt.tight_layout();
plt.savefig('conveergence-abs-w-m.pdf', format='pdf');
plt.show()
exit()







if 0:
	header = headers[i];
	f=open(fname,'ab');
	np.savetxt(f, a,fmt='%15.10f %15.10f %15.10f %15.10f', delimiter=' ',header=header,comments=' ')
	f.close();
	f=open(fname,'at')
	print('    ',file=f)
	print('    ',file=f)
	f.close();
	# waterfall plot ,comments='#'
	iy = i
	print('i == 1 ************')
	if show:
		offset = -iy*1 # (ny-iy)*5
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
	plt.ylim(-1*ignu-0,1);
	plt.xlabel("$\omega-\omega_0$",color='k',fontsize=24);
	plt.ylabel("$\Im G^R$",color='k',fontsize=24);
	plt.title(ttl,color='k',fontsize=16);
	plt.yticks([], [])
	plt.tight_layout();
	plt.savefig('abs-from-td-waterfall.pdf', format='pdf');
	plt.show();


