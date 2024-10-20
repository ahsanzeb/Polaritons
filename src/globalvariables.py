import numpy as np

#----------------------------------
inpfile = 0;
#----------------------------------

#----------------------------------
# set/calculated in readinput:
#----------------------------------
show = 0; corrtd= 0; matelem = 0;
groundstate = 0; justenergy = 0;
photonfraction=0;
zeroTPL = 0;
temp = 0;

diffoutdir = 0;
usenumpy = 0;
gamma=0; kappa=0;
e1=0; e2=0;
dt=0;
tf=0;
nwmax = 0; ntmax = 0;
printstep = 0;
nstates=0;

lamb0=0;
wr=0;wx=0;wc=0;wv=0;
n=0;m=0;mx=0;

Np=0;
lmin=0;lmax=0; nlmax=0;
loopover=0; 
lambda0=0;
lamlist=[];
lambin0 = [];
uselamlist = 0;

wclist = [];
wvlist = [];

itermax = 0;
tolr=0;

detuning=0;
eshft = 0;

ld=0;
dumy=0;

time0 = 0;
time00 = 0;
mem0 = 0;

#----------------------------------

#----------------------------------
# set/calculated in functions.countbasis()
n1fsym,n1,n2,n3,ntot=0,0,0,0,0;
listn1fsym=[];listn2 = []; listn3 = [];
#----------------------------------

#----------------------------------
# set/calculated in basis.fbasis():
MapNorms=np.zeros((1,1)); 
nlarge=False;
ndummy=5;

Nv1l=[];  
Norm1l=[]; Norm2l=[]; Norm3l=[];
map21=[]; 
Fact1l = [];
res = []; 
nrest = 0;
# set/used in mapping:
eloadlmap = [];
nblksizemap = [];
argxmap = [];
strtindsmap = [];
colsmap = [];
#----------------------------------

#----------------------------------
# set/calculated in hamiltonian:
Hcsm=[]; Hxsm=[]; Hvsm=[]; Hbsm=[]; Hgsm=[]; sft=[]; iden=[];
Hv0 = []; Hg=[];
#----------------------------------

#----------------------------------
# set/calculated in diagonalisation:
eigvv = [];
#----------------------------------

#----------------------------------
# matrixelements: defined global so that a function inside that module finds it! how weird!
evecs = [];
#----------------------------------

#----------------------------------
# multiple places: (matrixelements, gssolver, etc) to be used in eigensolver
#----------------------------------
ham1 = [];
lambin = []; lambin0 = []; 
ev0 = [];
#----------------------------------

nlist = [];
mlist = [];
psi0 = [];



