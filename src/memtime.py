
import globalvariables as o
import numpy as np

import time, psutil, os, sys

pid=os.getpid();
prno = psutil.Process(pid);
o.time00=time.time();
o.time0= 0;
# -------------------------------------------------
# Print time and memory information
def memtime(thismod):
	# get current time and memory usage:
	timnw=time.time() - o.time00;
	memnw = prno.memory_info().rss;
	# get the use between last and this call:
	tim= np.around(timnw-o.time0,2);
	mem = np.around((memnw-o.mem0)*1e-6,3);
	# print the message to stdout:
	x = ' total memory(MB)/time(s) = ';
	x += str(np.around(memnw*1e-6,3))+"/"+str(np.around(timnw,3));
	#with open('tout/'+thismod+'.txt','ab') as f_handle:
	#	np.savetxt(f_handle,[tim],fmt='%15.8f')
	if thismod != 'everything':
		x += ' used by '+thismod+' = '+str(mem)+"/"+str(tim);
	print(x)
	sys.stdout.flush()

	f=open('time-and-mem.txt','at')
	x = str(o.n)+' '+str(o.m)+' '+str(o.mx)+' ';
	x+= str(np.around(memnw*1e-6,3))+' '+str(np.around(timnw,3))+' ';
	x+= str(mem)+' '+str(tim);
	print(x,file=f)
	f.close()


	# reset the global var for later use:
	o.time0 = timnw; 	o.mem0 = memnw;
	return
# -------------------------------------------------


