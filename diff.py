import numpy as np

##### Diff mellan mina o antons resultat, datan ska va sparad som x column xp column
antonsdata = raw_input('Enter anton\'s datafile name:')
if len(antonsdata) < 1 : antonsdata = "100outanton.txt" # if no filename given, this file will be used
oscarsdata = raw_input('Enter oscar\'s datafile name:')
if len(oscarsdata) < 1 : oscarsdata = "out100oscarh.txt" # if no filename given, this file will be used

antonx, antonxp, antony, antonyp = np.loadtxt(antonsdata,unpack = True)
oscarx, oscarxp, oscary, oscaryp = np.loadtxt(oscarsdata,unpack = True)

diffx = antonx - oscarx
diffxp = antonxp - oscarxp
diffy = antony - oscary
diffyp = antonyp - oscaryp

stdx = np.std(diffx)
stdxp = np.std(diffxp)
stdy = np.std(diffy)
stdyp = np.std(diffyp)

print 'stdx:',stdx
print 'stdxp:',stdxp
print 'stdy:',stdy
print 'stdyp:',stdyp

# std for xin,xpin,yin,ypin
#print 'Initial beam std (when antonsdata is the init while and not results...)'
#print 'stdx:',np.std(antonx)
#print 'stdxp:',np.std(antonxp)
#print 'stdy:',np.std(antony)
#print 'stdyp:',np.std(antonyp)

## TODO: 
#1: make the program work by calling something like: python diff.py out.txt out2.txt