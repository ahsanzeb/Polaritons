
import numpy as np
import matplotlib.pyplot as plt


a=np.array([ 0.32751779,0.28948756,0.20891965,0.13057478])
b=-np.array([ -0.32751779,-0.40939724, -0.29545701, -0.18466063])
norm =np.array([1.0, 2.0, 2.0, 2.0])
c = range(4)

d = a*np.sqrt(norm);

print(d)

plt.plot(c,a,'sr',ms=10)
plt.plot(c,b,'sg',ms=10)
plt.plot(c,d,'ob',ms=5)

plt.show()

