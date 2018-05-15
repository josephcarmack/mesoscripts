#!/usr/bin/python
import paraview.simple as ps
import numpy as np
import paraview as pv
from vtk.util.numpy_support import vtk_to_numpy
import os
from sys import argv
import fnmatch
import curvature_analysis as curv

# main script

data_directory = argv[1]
save_directory = argv[2]
plot_name = argv[3]

filelist = os.listdir(data_directory)
filelist = fnmatch.filter(filelist,'fluid*')
filelist.sort(key=lambda x: int(x.split('.')[0].split('_')[1]))

# create ads array
points=len(filelist)
gc = np.zeros(points)
mc = np.zeros(points)
gc_per = np.zeros(points)
mc_per = np.zeros(points)
i = 0

# calculate gaussian and mean curvatures for each data file
for filename in filelist:
	a,b,c,d = curv.calc_curvature(data_directory+'/'+filename,'Viz',1.0,1.0,1)
	gc[i] = a
	mc[i] = b
	gc_per[i] = c
	mc_per[i] = d
	i = i + 1

# save data
skip = int(filelist[points-1].split('.')[0].split('_')[1])/(points-1)
t = np.arange(0,points*skip,skip) 

# save the ads data for replotting or manipulation later
np.savez(save_directory+'/'+plot_name+'.npz',gc,mc,gc_per,mc_per,t)
