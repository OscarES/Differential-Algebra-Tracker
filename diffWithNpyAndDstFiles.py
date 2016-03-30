import numpy as np

##### Diff between particle data, data should be stored as x column xp column y column yp
firstdata = raw_input('Enter first datafile name:')
if len(firstdata) < 1 : firstdata = "data/compareSC/" + "outpartFODSOspaceChargetest1.txt" # if no filename given, this file will be used
seconddata = raw_input('Enter second datafile name:')
if len(seconddata) < 1 : seconddata = "data/compareSC/" + "outpartFODSOspaceChargetest2.txt" # if no filename given, this file will be used

firstx, firstxp, firsty, firstyp = np.loadtxt(firstdata,unpack = True)
secondx, secondxp, secondy, secondyp = np.loadtxt(seconddata,unpack = True)

diffx = firstx - secondx
diffxp = firstxp - secondxp
diffy = firsty - secondy
diffyp = firstyp - secondyp

stdx = np.std(diffx)
stdxp = np.std(diffxp)
stdy = np.std(diffy)
stdyp = np.std(diffyp)

print 'stdx:',stdx
print 'stdxp:',stdxp
print 'stdy:',stdy
print 'stdyp:',stdyp

# std for xin,xpin,yin,ypin
#print 'Initial beam std (when firstsdata is the init while and not results...)'
#print 'stdx:',np.std(firstx)
#print 'stdxp:',np.std(firstxp)
#print 'stdy:',np.std(firsty)
#print 'stdyp:',np.std(firstyp)

## TODO: 
#1: make the program work by calling something like: python diff.py out.txt out2.txt