import numpy as np
import matplotlib.pyplot  as plt
from numpy import loadtxt
import os, fnmatch
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]




dos   = loadtxt("outdata/final/dos.dat", comments="#", unpack=False)
e     = loadtxt("outdata/final/E.dat", comments="#", unpack=False)
m     = loadtxt("outdata/final/M.dat", comments="#", unpack=False)
T     = loadtxt("outdata/final/T.dat", comments="#", unpack=False)
c     = loadtxt("outdata/final/c.dat", comments="#", unpack=False)
c_err = loadtxt("outdata/final/c_err.dat", comments="#", unpack=False)
d     = loadtxt("outdata/final/D.dat", comments="#", unpack=False)
d_err = loadtxt("outdata/final/D_err.dat", comments="#", unpack=False)
u     = loadtxt("outdata/final/u.dat", comments="#", unpack=False)
u_err = loadtxt("outdata/final/u_err.dat", comments="#", unpack=False)

cols = 16    #number of mpi threads
rows = 50000 #number of samples

# Use this to read all the energy-files under outdata/samples/*

# E = np.zeros((rows,cols))
# for cpu in range(0,cols):
#     directory = "outdata/samples/" + str(cpu) + "/"
#     files =  fnmatch.filter(os.listdir(directory),"energy_*")
#     # files.sort(key=natural_keys)
#     j = 0
#     for file in files:
#         if (j >= rows):
#             break
#         else:
#             energy = loadtxt(directory + file,unpack=False)
#             E[j][cpu] = energy
#             j = j+1
#

# Use this to store the energy values as a matrix to a single file
# with open('energies.txt','wb') as f:
#     np.savetxt(f, E, fmt='%d')

# Load the energies if you've used the code above previously
E = loadtxt('energies.txt', unpack=False)


plt.figure(0)
plt.plot(e,dos)
plt.figure(1)
plt.plot(T,c)
plt.figure(2)
plt.plot(T,u)
plt.figure(3)
plt.hist(E, bins=e, rwidth=0.8,normed=False)
plt.figure(4)
plt.hist(E.flatten(), bins=e, rwidth=0.8,normed=False)

plt.show()
