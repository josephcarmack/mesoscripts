#!/usr/bin/python

"""
This script takes the output of the bijel_ads.py script and fits a curve to it.
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from scipy.optimize import curve_fit

# curve fitting
def fitfunct(t,A,n,B):
    return A*t**n + B;
fitPar,fitCov = curve_fit(fitfunct,t,ads)
tt = np.linspace(skip,points*skip,num=1000)
lfit = fitPar[0]*tt**fitPar[1] + fitPar[2]
par_string = '$A =$ '+str(fitPar[0])+', $n =$ '+str(fitPar[1])+' B = '+str(fitPar[2])

# plotting
f = plt.figure(1)
ax = f.add_subplot(111)
plt.plot(t,ads,'ko',tt,lfit,'b')
plt.title('average domain size data with curve of the\nform $A*t^n$ fit to the data')
plt.xlabel('time in reduced units')
plt.ylabel('$l(t)$')
ax.annotate(par_string,xy=(t[int(t.size/2.0)],ads[int(ads.size)]/3.0))
plt.savefig(save_directory+'/'+plot_name+'.png')
