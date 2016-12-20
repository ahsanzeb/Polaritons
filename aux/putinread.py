

from scipy.sparse.linalg import lobpcg

import psutil
import gc

## from scipy import stats
import scipy
from scipy import stats
from scipy.sparse import identity
import numpy as np
import math
import os, sys
import time
## file basis.py in this folder
import basis
## import functions as fn
from multiprocessing import Pool
from scipy.integrate import ode
from scipy.fftpack import fft
import matplotlib.pyplot as plt
from numpy import sin, pi, exp
from scipy import arange
from scipy.sparse import coo_matrix



from readinput import * # get input variables from param.py or set default
import basis
from hamiltonian import * # fHb,fHbP2ss,fHbP2os,fHbP2osn2,fHbP12ss,fHb1,fHb2,fHb2n2,getHg,fHv1c,fHv0c,fHc,fHx
from diagonalisation import fdiagl, fdiagw;
from densitymatrices import * #getrho0,getrho1,getrho2,getrho2n2
from absorption import * #getCos2th, abs_matelem,phot_fraction,getwlist, fInteg, gettd, getFourierTransform, getFT,getFrewWindowCoor,fwriteCorr,fwriteFT,getFrankCondonEtc,Green
from functions import * # symmetrize, flat, mkchunks,


## start time
time0=time.time()

## for memory
pid=os.getpid();
prno = psutil.Process(pid);

#-------------------------------------
# to supress printing of small floats, print them 0
np.set_printoptions(suppress=True)
#-------------------------------------










# -------------------------------------------------
# Print time and memory information
timehamil=time.time()
tim= timehamil-timebasis
mem = prno.memory_info().rss;
print('      2. memory used now = '+str(mem*1e-6)+' MB   time taken = '+str(tim)+' s')
sys.stdout.flush()
# -------------------------------------------------


# -------------------------------------------------
# Print time and memory information
timebasis=time.time()
tim= timebasis - time0
mem = prno.memory_info().rss;
print('      1. memory used now  = '+str(mem*1e-6)+' MB   time taken = '+str(tim)+' s')
sys.stdout.flush()
t0old=timebasis
# -------------------------------------------------




ttot = time.time() - time0
print("   Total time taken = ",ttot,'s')
#****************************************************************

# for checking time taken by the prog vs N
f=open('timings-vs-n','ab')
np.savetxt(f, np.array([n,ttot,m,n1,n2,ntot])[None], fmt='%5i %15.10f %5i %5i %5i %5i')
f.close()

