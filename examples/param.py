
# ******************************************************
# for absorption 
# ******************************************************
#absorption = 'true';
#td = 'true';
#nstates = 10; # -ve means all states
#n=2;# n = no. of sites
#m=50;

nlist = [3]*2; # range(2,3,1); # [2,3,4,5]
# mlist=[[5,5],[5,40],[40,40]];
mlist=[[20,20],[20,25]]
ld =0.5; # amount in units of lamb0 to displace the basis

#mlist=[[5,5]]; # [[50,50],[30,40],[10,15],[6,10],[6,10],[6,10]];
#m=2; # cutoff on vib basis
#mx=3; #cuttoff for extra basis

if 0:
	mlist=[[5,5]];
	for m in range(1,100,1):
		mlist.append([5,5+m])

diffoutdir = 0; # old vsn='true';

# energies and couplings
wx=2; # exciton energy
wc=2; # cavity freq
wv=0.2; # vibrational freq

kappa = 0.05; # cavity's total decay rates def = 0.05; & k_R= k_L = k/2
gamma = 0.05; # molecules' decay rates def = 0.05
tf = 100.0; # total time tf in units of 1/wr; def = 100
dt = 0.1 # time interval for correlation data; TDSE is integrated between all t,t+dt
nwft = 300;

e1,e2 = -1.5,2.5; # energy window for FT absorption data, w.r.t. wx=0
#show = 1; # 0,1 show plot of absorption and correlation functions
# ******************************************************

# ******************************************************
# for ground state calculations
# ******************************************************
# stop after energy calculations. Dont calc dm's.
#onlyenergy = 'true'

# Now after loop over lambda/wr and N:
# wr and lamb0 are defined by the old style for loopover wr/lambda 
# ###wr=0.5; # rabi freq; wr = g/sqrt(n)
# ###lamb0 = 0; # lambda0

wr=0.5; # wr = g/sqrt(n), for loop over lambda0 calculations
wmin = 0.5
wmax = 2 
nwmax = 1 # for delta_wr = 0.25

# nwmax=nwft # old usage: FT no. of freq points

lambda0 = 1.5; # lambda0 for loop over wr calculations
lmin=2.0;
lmax=2.5;
nlmax=1 # for delta_wr = 0.25

# set loop over either lambda0 or wr:
loopover = 'lambda0' 	#  loop over lambda0 (\lambda_0)
#loopover = 'wr'  	#  loop over wr (\omega_R)

#lamb0 = lmin; # old lam0, just for compatibility, hope not needed anymore
# ******************************************************



# ******************************************************
# for any calculation
# ******************************************************
# ******************************************************
# for both absorption (if td != true) and gs calculations
# ******************************************************
Np = 8; # number of processors
# options for lobpcg
itermax = 200; 
tolr=1e-10;
# ******************************************************


