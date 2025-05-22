import numpy as np
import math
import scipy.integrate as integrate
import h5py
import os
from constants import *
from config import *


def getdomain(file):
    infile = open(file)
    lines = infile.readlines()
    dom_range = np.zeros((2,3))
    ncell = np.zeros(3)
    dom_min = [0.0,0.0,0.0]
    dom_min[0] = float(lines[3].split()[2])
    dom_min[1] = float(lines[3].split()[3])
    dom_min[2] = float(lines[3].split()[4])
    
    dom_max = [0.0,0.0,0.0]
    dom_max[0] = float(lines[4].split()[2])
    dom_max[1] = float(lines[4].split()[3])
    dom_max[2] = float(lines[4].split()[4])

    ncell[0]=int(lines[15].split()[2])
    ncell[1]=int(lines[15].split()[3])
    ncell[2]=int(lines[15].split()[4])
    
    return dom_min, dom_max, ncell


def get_folder_size(directory):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            # Skip if it is symbolic link
            if not os.path.islink(file_path):
                total_size += os.path.getsize(file_path)
    return total_size


def get_temp(rho_gas, eint):
    
    file = '/Users/aditivijayan/Projects/SKA-obs/data/CloudyData_UVB=FG2011.h5'
    grackle = h5py.File(file)
    array = grackle['CoolingRates/Primordial/MMW'][()]

    table = array[:,0,:]
    table_nH   = np.logspace(-10., 4, array.shape[0])
    table_temp = np.logspace(1,  9, array.shape[2])
    
   
    
    i=0
    bins = 200
    egas_arr = np.logspace(-25., -1., bins)
    nH_arr   = np.logspace(-7.0, 4.0, int(bins))
    T = np.zeros((egas_arr.shape[0],nH_arr.shape[0]))

    for egas in egas_arr:
        j=0
        for nH in nH_arr:
            C = (gamma - 1.) * egas / (boltzmann_constant_cgs*nH)
            minT = C*np.amin(table)
            maxT = C*np.amax(table)
            def func(T):
                mu = interpolate.interp2d(table_temp, table_nH, table,\
                                  kind='linear', copy=True, bounds_error=False, fill_value=None)
                return C*mu(T,nH)[0] - T

            T[i,j] = scipy.optimize.toms748(func, minT, maxT)
            j+=1
        i+=1

    egas0=eint
    density = rho_gas
    cloudy_H_mass_fraction = 1. / (1. + 0.1 * 3.971)
    rho0 = density*cloudy_H_mass_fraction/hydrogen_mass_cgs


    logrho_arr = np.log10(nH_arr[:-1])
    logrho     = np.log10(rho0)
    delta_rho  = logrho_arr[1] - logrho_arr[0]
    idxrho     = (np.floor((logrho - np.amin(logrho_arr))/delta_rho)).astype('int')

    logEgas_arr = np.log10(egas_arr[:-1])
    logEgas     = np.log10(egas0)
    delta_egas  = logEgas_arr[1] - logEgas_arr[0]

    idxegas     = (np.floor((logEgas-np.amin(logEgas_arr))/delta_egas)).astype('int')


    wgt_rho  = (logrho - (np.amin(logrho_arr) + delta_rho*idxrho))/delta_rho
    wgt_egas = (logEgas - (np.amin(logEgas_arr) + delta_egas*idxegas))/delta_egas

    temp = (1.-wgt_rho)*(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho)]   +\
               wgt_rho *    wgt_egas * T[tuple(idxegas+1), tuple(idxrho+1)] +\
          (1. -wgt_rho)*    wgt_egas * T[tuple(idxegas+1), tuple(idxrho)]   +\
               wgt_rho *(1.-wgt_egas)* T[tuple(idxegas)  , tuple(idxrho+1)]  
    
    return temp

def limits(arr):
    return np.amax(arr), np.amin(arr)
