# ------------------------------
# photon block:
def fHb1(jchunk):
	hb1 = [];
	for jj in range(jchunk[0],jchunk[1]):
		m1 = Norm2l[jj];
		for mj1 in range(0,m):
				mj2 = mj1 + 1;
				i1 = map21[jj,mj1];
				i2 = map21[jj,mj2];
				m11=Norm1l[i1];
				m12=Norm1l[i2];
				yy = math.sqrt(mj2);
				xx = math.sqrt(m11*m12);
				xij = -ld*m1*n*yy/xx; # -ld displaced 
				hb1.append([i1,i2,xij]);
	return hb1
# ------------------------------
# exciton block: excited-site-diagonal part
# ------------------------------
def fHb2(kchunk):
	hb2 = [];
	for kk in range(kchunk[0],kchunk[1]):
		m1 = Norm3l[kk];
		for mj1 in range(0,m):
				mj2 = mj1 +1;
				i1 = map32[kk,mj1];
				i2 = map32[kk,mj2];
				m11=Norm2l[i1];
				m12=Norm2l[i2];
				yy = math.sqrt(mj2);
				xx = math.sqrt(m11*m12);
				xij=-ld*m1*(n-1)*yy/xx; # -ld displaced
				for m2 in range(0,m+1):
					ii1=m2*n2 + i1; # diag in m2 contrib only
					ii2=m2*n2 + i2;
					hb2.append([n1+ii1,n1+ii2,xij]);
	return hb2
# ------------------------------	
def fHb2n2(mj1):
	hb2 = [];
	mj2 = mj1+1;
	i1 = mj1;
	i2 = mj2;
	yy = math.sqrt(mj2);
	xij=-ld*yy; # -ld displaced
	for m2 in range(0,m+1):
		ii1=m2*n2 + i1; # diag in m2 contrib only
		ii2=m2*n2 + i2;
		hb2.append([n1+ii1,n1+ii2,xij]);
	return hb2
# ------------------------------	

