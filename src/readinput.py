

import globalvariables as o
import os, sys
import numpy as np
from numpy import pi

#888888888888888888888888888
# system parameters
#888888888888888888888888888

# path to input parameter file
path_param = os.getcwd();

# add this path to python search path
sys.path.insert(0, path_param)

# import param
# import the input parameters
inpfname = o.inpfile;
inpfname = os.path.splitext(o.inpfile)[0]
import importlib
param = importlib.import_module(inpfname);
#-------------------------------------
# read input file and set default values
# search input file for paramters if they are needed.
#-------------------------------------
print(' ')# empty line

print(' ====> reading input file ... ')


show = 0; corrtd= 0; matelem = 0;
groundstate = 0;
justenergy = 0;
photonfraction=0;
nstates=1;
absorption='false';
td = 'false';
wr = 0; lambda0 = 0; 
lamlist=[]; uselamlist = 0; 
lambin0 = [];
lmin = 0;
lmax = 0;
nlmax = 1;
# loopover = 'lambda0';
wclist = [];  # wclist, wvlist etc have an issue... works only for a single value of n
# (readinout: len(nlist)==len(wclist), works otherwise default wc=2 is used!!! )
# (this means we cant use wclist for a nlist with multiple n values... 
#Solution: change the style of loopover... include all loops there... 
#over wc, wv,  etc... and use a single function for diagonalisation
# instead of currect style two seperate for lambda and wr.)
wvlist = [];

temp = 0;
if hasattr(param, 'temp'):
	temp = param.temp
	#print(' Finite Temperature for absorption... ')
	

if hasattr(param, 'wclist'):
	if hasattr(param, 'nlist'):
		if len(param.wclist)==len(param.nlist):
			wclist = param.wclist
			print(' vs wc (detuning): wclist will be used... ')
if hasattr(param, 'wvlist'):
	if hasattr(param, 'nlist'):
		if len(param.wvlist)==len(param.nlist):
			wvlist = param.wvlist
			print(' vs wv: wvlist will be used... ')

if hasattr(param, 'nlarge'):
	o.nlarge=param.nlarge
	print('readinp: nlarge = ',o.nlarge)
	if(o.nlarge):
		print('readinp: Basis with permutations for N='+str(o.ndummy)+' will be used.... !')

			
if hasattr(param, 'absorption'):
	absorption = param.absorption
	if (absorption=='true'):
		print(" absorption = TRUE ");
		# .............................................
		if hasattr(param, 'td'):
			td = param.td
			if (td=='true'):
				print(" td = TRUE ");				
				nstates = 1;
				if hasattr(param, 'gamma'):
					gamma = param.gamma;
				else:
					gamma = 0.05;
					print(" gamma not set: default = 0.05 will be used...");
				if hasattr(param, 'kappa'):
					kappa = param.kappa;
				else:
					kappa = 0.05;
					print(" kappa not set: default = 0.05 (& k_L = k/2) will be used...");
				if hasattr(param, 'e1') and hasattr(param, 'e2') :
					e1 = param.e1; e2 = param.e2;
				else:
					e1,e2 = -2,4;
					print(" e1,e2 not set: default = -2,4 will be used...");
				if hasattr(param, 'dt'):
					dt= param.dt;
					if dt < 0.005:
						dt = 0.01;
						print(" dt given < 0.005: setting to 0.01 ...");
				else:
					dt = 0.01;
					print(" dt not set: default = 0.01 will be used...");
				if hasattr(param, 'tf'):
					tf= param.tf;
				else:
					tf = 100;
					print(" tf not set: default = 100 will be used...");
					print("    ==> ntmax ~ tf/dt = ", int(tf/dt) );
				# --------------------------
				tdecay = 30/(kappa+gamma); # time corr decay to ~ 0
				if tf > tdecay:
					#tf = tdecay;
					#print(" tf > 30/(kappa+gamma): setting tf = 30/(kappa+gamma) = ",tdecay);
					print(" tf > 30/(kappa+gamma):  tf, 1/(kappa+gamma)= ",tf,1/(kappa+gamma));
					print(" tf ~ 30/(kappa+gamma) should be ok ... ");
				# Fourier Transform: 
				# freq range within one 'period'.
				Wft = 2*pi/dt; # period on freq axis
				if e2-e1 > Wft:
					e1 = -Wft/2;
					e2 = +Wft/2;
					print(" e2-e1 > 2*pi/dt: setting e1=-pi/dt, e2=pi/dt = ",e1,e2);
				# --------------------------
				if hasattr(param, 'nwft'):
					nwmax = param.nwft;
				else:
					nwmax = 200;
					print(" nwft not set: default = 200 will be used...");
				if hasattr(param, 'show'):
					show= param.show;
			else:
				td = None;
		else:
			td = None;
		# .............................................
		if td !='true':
			print("NOT td: eigenstates matrix elemets will be calculated ...");			
			if hasattr(param, 'nstates'):
				nstates= param.nstates;
			else:
				nstates = 400;
				print(" nstates not set: default = 400 will be used...");
	# .............................................
	else:
		absorption = None
else:
	absorption=None
# .............................................

printstep = 25;
if hasattr(param, 'printstep'):
	printstep = param.printstep

usenumpy = 0;
if hasattr(param, 'usenumpy'):
	usenumpy = param.usenumpy
	if usenumpy != 0 and usenumpy !=1:
		print(' set usenumpy=0/1 for numpy/fortran for time evolution')
		exit()
o.usenumpy = usenumpy;

zeroTPL = 0;
if hasattr(param, 'zeroTPL'):
	zeroTPL = param.zeroTPL
	if zeroTPL != 1 and zeroTPL !=0:
		print(' Absorption or PL? zeroTPL !=0,1; set it properly...' )
o.zeroTPL = zeroTPL


if hasattr(param, 'photonfraction'):
	photonfraction = param.photonfraction

#------------------------------------- 
if hasattr(param, 'onlyenergy'):
	onlyenergy = param.onlyenergy
	if (onlyenergy=='true' and absorption=='true'):
		onlyenergy = None;
	elif(onlyenergy=='true' and absorption !='true'):
		print("       onlyenergy = TRUE ")
	else:
		onlyenergy = None
else:
	onlyenergy=None
#-------------------------------------
if(absorption=='true' and onlyenergy=='true'):
	onlyenergy=None;
	#print("       absorption=T so setting onlyenergy=F")
#-------------------------------------
diffoutdir = 1;
if hasattr(param, 'diffoutdir'):
	diffoutdir = param.diffoutdir
elif absorption=='true':
	diffoutdir =0;
	print(" For absorption, default diffoutdir=0 will be used...");


print(" ")
print("Warning: For DM calculations, mx=m only implemented! ");
print(" ")

#-------------------------------------
# set some default parameters:
if 1:
	if hasattr(param, 'wx'):
		wx= param.wx;
	else:
		wx = 2.0;
		print(" wx not set: default = 2.0 will be used...");
	if hasattr(param, 'wc'):
		wc= param.wc;
	else:
		wc = 2.0;
		print(" wc not set: default = 2.0 will be used...");
	if hasattr(param, 'wv'):
		wv= param.wv;
	else:
		wv = 0.2;
		print(" wv not set: default = 0.2 will be used...");

	if hasattr(param, 'nlist'):
		nlist = param.nlist;
		n = nlist[-1]; # last, probably the largest
	elif hasattr(param, 'n'):
		n= param.n;
	else:
		n = 1; nlist = [1];
		print(" n not set: default = 1 will be used...");	

	# if lamlist, then use for m in mlist, lambda0 in lamlist,
	if hasattr(param, 'lamlist'):
		lamlist = param.lamlist;
		if len(nlist)==1:
			nlist = nlist*len(lamlist);

	if hasattr(param, 'mlist'):
		mlist = param.mlist;
		if isinstance(mlist[0], list):
			mlast = mlist[-1];
		else:
			mlast = mlist; # last, probably the smallest
			mlist=[mlist];
	elif n<4:
			mlast =[10,10]; mlist=[mlast];
			print(" mlist not set: default [10,10] (n<4) will be used...");
	else:
			mlast =[2,2]; mlist=[mlast];
			print(" mlist not set: default [2,2] (n>4) will be used...");
	if len(mlist) < len(nlist):
		for i in range(len(nlist)-len(mlist)):
			mlist.append(mlast);

	if hasattr(param, 'Np'):
		Np= param.Np;
	else:
		Np = 8;
		print(" Np not set: default = 8 will be used...");


if hasattr(param, 'loopover'):
	loopover = param.loopover;
	# loopover which paramter, lambda0 or wr?
	if loopover == 'lambda0':
		if hasattr(param, 'lamlist'):
			uselamlist = 1;
			lamlist = param.lamlist;
			nlmax = 1; # do one by one, as if a nlist is given with nlmax=1.
			print('   lamlist is given... warning: lambda0 loop not parallel!')
		else:
			uselamlist = 0;
			wr=param.wr # g/sqrt(n)
			lmin=param.lmin;
			lmax=param.lmax;
			nlmax=param.nlmax;
			lambin0 = np.linspace(lmin,lmax, nlmax);
			print('   loopover == "lambda0" ... ')
	elif loopover == 'wr':
		lambda0 = param.lambda0
		lmin=param.wmin;
		lmax=param.wmax;
		nlmax=param.nwmax;
		lambin0 = np.linspace(lmin,lmax, nlmax);
		print('   loopover = "wr" ... ')
	else: # if absorption !='true':
		print('   error: set option loopover to lambda0 or wr (with relevant min,max,nmax data)')
		exit()







	
#-------------------------------------
#  set some default parameters for lobpcg
if td != 'true':
	if hasattr(param, 'itermax'):
		itermax= param.itermax;
	else:
		itermax = 150;
		print(" itermax not set: default = 150 will be used...");
	if hasattr(param, 'tolr'):
		tolr= param.tolr;
	else:
		tolr = 1e-8;
		print(" tolr not set: default = 1e-8 will be used...");
#-------------------------------------
if (abs(wc-wx)<1e-6 and td != 'true'):
		detuning = 0;
		# print("			detuning < 1e-6")
		eshft= wx; # if eshft < 1e-6 then consider it 
		wx=0; wc=0;
else:
		detuning = 1;
		eshft = 0;

#-------------------------------------
# vib basis are displaced in {x_i} by ld*lamb0
if hasattr(param, 'ld'):
	ld = param.ld;
	lamtol = 0.5; 
	if (loopover=='wr' and lambda0<lamtol): ld = 0;
	if loopover=='lambda0':
		if (uselamlist and len(lamlist)==1 and lamlist[0]<lamtol): ld = 0;
		elif ( not uselamlist and nlmax==1 and lambin0[0]<lamtol): ld = 0;
	if (ld > 0.5 or ld < 0):
		print(" only 0 < ld < 0.5 can be useful. (no? ld>0.5 better if state is excitonic!)")
		#exit()
else:
	ld=0.5;
#-------------------------------------

#****************************************************************
# set global variables:
#****************************************************************

# which type of calculation?
if absorption=='true':
	if td == 'true':
		corrtd = 1;
	else: 
		matelem = 1;  o.matelem = 1;
else:
	groundstate = 1;
	if onlyenergy == 'true':
		justenergy = 1;
	#else:
	#	print(' vib dm: setting mx=m, mx > m not implemented... ');
	#	mx = m;

# for all types of calculations
o.nlist, o.mlist = nlist, mlist;
# o.n, o.m, o.mx, o.Np =n,m,mx,Np;
o.Np = Np;

o.wclist = wclist;
o.wvlist = wvlist;

o.wr, o.wx, o.wc, o.wv=wr,wx,wc,wv
o.ld=ld
o.nstates = nstates
o.diffoutdir = diffoutdir

# cheating: apply loopover lambda0 to Ground vib band instead
#n=min(nlist);m=min(mlist)[0]
#if n>1 and m>2 and temp ==1 and loopover=='lambda0':
if temp ==1 and loopover=='lambda0':
	o.temp= temp;
	lmin,lmax, nlmax = lmin,lmin,6; # 6 initial states needed...
	lambin0 = np.linspace(lmin,lmax, nlmax);
	print(' temp=1 so cheating.... lmin will be used...')
	print(' Finite Temperature for absorption... for 6 states with Nv <= 3 as the starting states for time evolution')
	
o.lmin, o.lmax,  o.nlmax = lmin,lmax, nlmax
o.loopover =loopover; 
o.lambda0=lambda0
o.lambin0 = lambin0;
o.uselamlist = uselamlist;
o.lamlist = lamlist;

# for specific type of calculations
if corrtd:
	o.corrtd = corrtd;
	o.gamma, o.kappa = gamma, kappa
	o.e1, o.e2, o.dt, o.tf = e1,e2,dt,tf
	o.nwmax = nwmax ;
	o.show = show
	o.ntmax = int(tf/dt);
	o.printstep = printstep;
	
if matelem or groundstate:
	o.itermax,  o.tolr=itermax, tolr
	o.detuning, o.eshft = detuning, eshft

if groundstate:
	o.groundstate = groundstate;
	o.photonfraction = photonfraction;
	if justenergy:
		o.justenergy = justenergy;
#****************************************************************


