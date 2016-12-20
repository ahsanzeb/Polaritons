
#----------------------------------
inpfile = 0;
#----------------------------------

#----------------------------------
# set/calculated in readinput:
#----------------------------------
show = 0; corrtd= 0; matelem = 0;
groundstate = 0; justenergy = 0;

gamma=0; kappa=0;
e1=0; e2=0;
dt=0;
tf=0;
nwmax = 0;
nstates=0;

lamb0=0;
wr=0;wx=0;wc=0;wv=0;
n=0;m=0;mx=0;
Np=0;
mg=0;
lmin=0;lmax=0; nlmax=0;
loopover=0; 
lambda0=0;

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
listn2 = []; listn3 = [];
#----------------------------------

#----------------------------------
# set/calculated in basis.fbasis():
Nv1l=[]; Nv2l=[]; 
Norm1l=[]; Norm2l=[]; Norm3l=[];
map21=[]; map32=[];
#----------------------------------

#----------------------------------
# set/calculated in hamiltonian:
Hcsm=[]; Hxsm=[]; Hvsm=[]; Hbsm=[]; Hgsm=[]; sft=[]; iden=[];
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


