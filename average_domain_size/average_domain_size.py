import numpy as np
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
import os
from sys import argv
import fnmatch

######################################################################
# functions
######################################################################

def calc_ads(filename):

    # read file
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(filename)
    reader.Update()
    data = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(0))

    # manipulate data
    shape = reader.GetOutput().GetDimensions()
    data = data.reshape(shape)
    spacing = reader.GetOutput().GetSpacing()
    dx,dy,dz = spacing
    nx,ny,nz = data.shape
    nxyz = nx*ny*nz

    # build the k-matrix
    kx1 = np.mod( 1.0/2 + np.arange(float(nx))/nx,1.0) -1.0/2
    kx = kx1*(2.0*np.pi/dx)
    ky1 = np.mod( 1.0/2 + np.arange(float(ny))/ny,1) -1.0/2
    ky = ky1*(2.0*np.pi/dy)
    kz1 = np.mod( 1.0/2 + np.arange(float(nz))/nz,1) -1.0/2
    kz = kz1*(2.0*np.pi/dz)
    KX,KY,KZ = np.meshgrid(kx,ky,kz)

    # calculate the average domain size
    phi_fft = np.fft.fft2(data)/nxyz
    Sk = np.abs(phi_fft)*np.abs(phi_fft)
    kx2 = np.sum(np.sum(KX*KX*Sk))/np.sum(np.sum(Sk))
    ky2 = np.sum(np.sum(KY*KY*Sk))/np.sum(np.sum(Sk))
    kz2 = np.sum(np.sum(KZ*KZ*Sk))/np.sum(np.sum(Sk))
    lx = 2.0*np.pi/kx2**0.5
    ly = 2.0*np.pi/ky2**0.5
    lz = 2.0*np.pi/kz2**0.5
    ads = (lx + ly + lz)/3.0
    return ads;
