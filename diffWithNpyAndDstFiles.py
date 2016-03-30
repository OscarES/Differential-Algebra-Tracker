import numpy as np
from IOHandler import loadMultipart

##### Diff between particle data, data should be stored as x column xp column y column yp
firstdata = raw_input('Enter first multipart datafile name:')
seconddata = raw_input('Enter second multipart datafile name:')

firstmultipart = loadMultipart(firstdata)
secondmultipart = loadMultipart(seconddata)

xf = [firstmultipart[i][0][0] for i in xrange(len(firstmultipart))]
xpf = [firstmultipart[i][0][1] for i in xrange(len(firstmultipart))]
yf = [firstmultipart[i][0][2] for i in xrange(len(firstmultipart))]
ypf = [firstmultipart[i][0][3] for i in xrange(len(firstmultipart))]
zf = [firstmultipart[i][0][4] for i in xrange(len(firstmultipart))]
zpf = [firstmultipart[i][0][5] for i in xrange(len(firstmultipart))]

xs = [secondmultipart[i][0][0] for i in xrange(len(secondmultipart))]
xps = [secondmultipart[i][0][1] for i in xrange(len(secondmultipart))]
ys = [secondmultipart[i][0][2] for i in xrange(len(secondmultipart))]
yps = [secondmultipart[i][0][3] for i in xrange(len(secondmultipart))]
zs = [secondmultipart[i][0][4] for i in xrange(len(secondmultipart))]
zps = [secondmultipart[i][0][5] for i in xrange(len(secondmultipart))]

diffx = np.array(xf) - np.array(xs)
diffxp = np.array(xpf) - np.array(xps)
diffy = np.array(yf) - np.array(ys)
diffyp = np.array(ypf) - np.array(yps)
diffz = np.array(zf) - np.array(zs)
diffzp = np.array(zpf) - np.array(zps)

diffx = diffx.astype('float')
diffxp = diffxp.astype('float')
diffy = diffy.astype('float')
diffyp = diffyp.astype('float')
diffz = diffz.astype('float')
diffzp = diffzp.astype('float')

stdx = np.std(diffx)
stdxp = np.std(diffxp)
stdy = np.std(diffy)
stdyp = np.std(diffyp)
stdz = np.std(diffz)
stdzp = np.std(diffzp)

print 'stdx:',stdx
print 'stdxp:',stdxp
print 'stdy:',stdy
print 'stdyp:',stdyp
print 'stdz:',stdz
print 'stdzp:',stdzp

# std for xin,xpin,yin,ypin
#print 'Initial beam std (when firstsdata is the init while and not results...)'
#print 'stdx:',np.std(firstx)
#print 'stdxp:',np.std(firstxp)
#print 'stdy:',np.std(firsty)
#print 'stdyp:',np.std(firstyp)

## TODO: 
#1: make the program work by calling something like: python diff.py out.txt out2.txt