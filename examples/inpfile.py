
# ******************************************************
# for absorption 
# ******************************************************
# absorption = 'true';
# td = 'true';
# nstates = 400; # -ve means all states
n=2;# n = no. of sites
nlist = [2,3,4,5,6,7];

m=2; # cutoff on vib basis
mx=4; #cuttoff for extra basis

# energies and couplings
wx=2; # exciton energy
wc=2; # cavity freq
wv=0.2; # vibrational freq

kappa = 0.05; # cavity's total decay rates def = 0.05; & k_R= k_L = k/2
gamma = 0.05; # molecules' decay rates def = 0.05
tf = 10.0; # total time tf in units of 1/wr; def = 100
dt = 1.0 # time interval for correlation data; TDSE is integrated between all t,t+dt
nwft = 30;
e1,e2 = -4,4; # energy window for FT absorption data, w.r.t. wx=0
#show = 1; # 0,1 show plot of absorption and correlation functions
# ******************************************************

# ******************************************************
# for ground state calculations
# ******************************************************
# stop after energy calculations. Dont calc dm's.
# onlyenergy = 'true'

# Now after loop over lambda/wr and N:
# wr and lamb0 are defined by the old style for loopover wr/lambda 
# ###wr=0.5; # rabi freq; wr = g/sqrt(n)
# ###lamb0 = 0; # lambda0

wr=1; # wr = g/sqrt(n), for loop over lambda0 calculations
wmin = 0
wmax = 3 
nwmax = 6 # for delta_wr = 0.25

lambda0 = 2.5; # lambda0 for loop over wr calculations
lmin=0;
lmax=1;
nlmax=2 # for delta_wr = 0.25

# set loop over either lambda0 or wr:
loopover = 'lambda0' 	#  loop over lambda0 (\lambda_0)
#loopover = 'wr'  	#  loop over wr (\omega_R)

# ******************************************************



# ******************************************************
# for any calculation
# ******************************************************
ld = 0.5; # amount in units of lamb0 to displace the basis
# ******************************************************
# for both absorption (if td != true) and gs calculations
# ******************************************************
Np = 8; # number of processors
# options for lobpcg
itermax = 200; 
tolr=1e-8;
# ******************************************************


