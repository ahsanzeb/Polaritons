

import globalvariables as o
import os, sys
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
nstates=1;
absorption='false';
td = 'false';
wr = 0; lambda0 = 0

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
				tdecay = 10/(kappa+gamma); # time corr decay to ~ 0
				if tf > tdecay:
					tf = tdecay;
					print(" tf > tdecay: setting tf = tdecay = 10/(kappa+gamma) = ",tdecay);
				if tf/dt > 10000:
					print(" to avoid ntmax > 10000: set tf/dt to default 100, 0.01 ...")
					if dt >= 0.01:
						tf = 1000*dt;
						print("ntmax > 10000; dt >= 0.01: setting tf = 1000*dt = ",tf)
					else:
						dt = tf/1000;
						print("ntmax > 10000; tf <= 100: setting dt = tf/1000 = ",dt)
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

#------------------------------------- 
if hasattr(param, 'onlyenergy'):
	onlyenergy = param.onlyenergy
	if (onlyenergy=='true'):
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
if hasattr(param, 'vsn'):
	vsn = param.vsn
	if (vsn=='true'):
		print("       vsn = TRUE ")
	else:
		vsn = None
		#print("       vsn = None ")
else:
	if absorption=='true':
		print(" For absorption, default vsn='true' will be used...");
		vsn='true';
	else:
		vsn=None
#-------------------------------------
# set some default parameters:
if 1:
	if hasattr(param, 'wr'):
		wr= param.wr;
	else:
		wr = 0.5;
		print(" wr not set: default = 0.5 will be used...");
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
	if hasattr(param, 'n'):
		n= param.n;
	else:
		n = 2;
		print(" n not set: default = 2 will be used...");
	if hasattr(param, 'm'):
		m= param.m;
	else:
		if n<4:
			m = 10;
		else:
			m = 2;
		print(" m not set: default = 10 (n<4) or 2 will be used...");
	if hasattr(param, 'mx'):
		mx= param.mx;
	else:
		if n<4 and m<5:
			mx = 10;
		else:
			mx = m;
		print(" mx not set: default = 10 (n<4,m<5) or m will be used...");
	if hasattr(param, 'Np'):
		Np= param.Np;
	else:
		Np = 8;
		print(" Np not set: default = 8 will be used...");



if(absorption=='true'):
	if hasattr(param, 'lamb0'):
		lamb0= param.lamb0;
	else:
		lamb0 = 1.0;
		print(" lamb0 not set: default = 1.0 will be used...");
	# set upper cutoff on vib spectrum for green function; ~ 10, not more. 
	# if mx > 10:
	# 	mg = 10;
	# else:
	# 	mg = mx;
	mg = 10;
	# set lmin etc to be consistent with eigensolver functions.
	lmin = lamb0;
	lmax = 2*lamb0;
	nlmax = 1;
	loopover = 'lambda0';


if hasattr(param, 'loopover'):
	loopover = param.loopover;
	# loopover which paramter, lambda0 or wr?
	if loopover == 'lambda0':
		wr=param.wr # g/sqrt(n)
		lmin=param.lmin;
		lmax=param.lmax;
		nlmax=param.nlmax;
		print('   calculations would be for loop over lambda0 ... ')
	elif loopover == 'wr':
		lambda0 = param.lambda0
		lmin=param.wmin;
		lmax=param.wmax;
		nlmax=param.nwmax;
		print('   calculations would be for loop over wr ... ')
	else:
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
print(" n, m, mx = ",n,m,mx);
if mx < m:
	mx=m; print("was mx < m =====> setting mx=m")
#-------------------------------------
if n==1:
	print('stop: n=1 not in python code..., use mathematica code')
	exit()
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
	if (absorption=='true' and lamb0 ==0):
		ld = 0;
else:
	ld=0
print(" ld = ", ld)
if (ld > 0.5 or ld < 0):
	print("  stop: only 0 < ld < 0.5 can be useful!")
	exit()

#-------------------------------------






#****************************************************************
# create dir for output dm
if vsn=='true':
	dumy="data/n-all"
else:
	dumy="data/n-"+str(n)
path = os.getcwd();
path = path+"/"+dumy
os.makedirs(path, exist_ok=True)




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
	else:
		print('vibrational dm: setting mx=m, mx > m not implemented... ');
		mx = m;


# for all types of calculations
o.n, o.m, o.mx, o.Np =n,m,mx,Np;
o.wr, o.wx, o.wc, o.wv=wr,wx,wc,wv
o.ld=ld
o.nstates = nstates
o.dumy=dumy

if 1:
	o.lmin, o.lmax,  o.nlmax = lmin,lmax, nlmax
	o.loopover =loopover; 
	o.lambda0=lambda0

# for specific type of calculations
if corrtd:
	o.corrtd = corrtd;
	o.gamma, o.kappa = gamma, kappa
	o.e1, o.e2, o.dt, o.tf = e1,e2,dt,tf
	o.nwmax = nwmax ; o.mg  = mg;
	o.show = show

if matelem or corrtd:
	o.lamb0 = lamb0;
	o.ntmax = int(tf/dt);


if matelem or groundstate:
	o.itermax,  o.tolr=itermax, tolr
	o.detuning, o.eshft = detuning, eshft

if groundstate:
	o.groundstate = groundstate;
	if justenergy:
		o.justenergy = justenergy;
#****************************************************************


