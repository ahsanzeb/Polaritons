



import globalvariables as o


#****************************************************************
# make even chunks
# divide range(0,nnn) in Npp (= # of processes) intervals with boundaries [i1,i2] in array listn
def mkchunks(Npp,nnn):
	n21 = int(nnn/Npp); rem2 = nnn % Npp;
	listn=[]; i1=0;
	for i in range(rem2):
		i2= i1 + n21 +1
		listn.append([i1,i2])
		i1=i2
	for i in range(rem2,Npp):
		i2= i1 + n21
		listn.append([i1,i2])
		i1=i2
	return listn
#****************************************************************

def countbasis():
	n,m,mx = 	o.n,o.m,o.mx;
	# counting basis states of various sectors:
	n1fsym=len(o.Norm1l); # old thing with name n1
	n2=len(o.Norm2l); # No. of symmetrised basis for |0c,1x>
	# new n1: so that we dont have change shifts by n1 everywhere
	n1 = n1fsym + (mx-m)*n2;
	# total number of basis states
	ntot = n1 + (mx+1)*n2;
	if (n>2):
		n3=len(o.Norm3l); # basis size for N-2 sites
		o.n1fsym,o.n1,o.n2,o.n3,o.ntot = n1fsym,n1,n2,n3,ntot; 
	o.n1fsym,o.n1,o.n2,o.ntot = n1fsym,n1,n2,ntot;
	return


def mkpoollists():
	n,n2,n3,Np = 	o.n,o.n2,o.n3,o.Np 
	o.listn2 = mkchunks(Np,n2);
	if (n>2):
		o.listn3 = mkchunks(Np,n3);
	return


