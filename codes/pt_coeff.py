
from classy import Class
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from bcoeff import bcoeff
from ls_coeff import lscoeff
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys




def ptcoeff(self, data, kclass, Massbins):
       
	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	####################################################################
	#### Store the mass bins available in an array
	mbins = ['M1', 'M2', 'M3', 'M4']

	
	####################################################################
	#### get the coefficients from the dat files in the data directory 
	#### get the large scale amplitude of bcc at different z and for different neutrino masses
	self.data_directory = data.path['root']

	Aprime = np.zeros((350,l2))
	Bprime = np.zeros((350,l2))
	Cprime = np.zeros((350,l2))
	Dprime = np.zeros((350,l2))
	Eprime = np.zeros((350,l2))
	Fprime = np.zeros((350,l2))
	Gprime = np.zeros((350,l2))
	Hprime = np.zeros((350,l2))
	Pmod_dd_prime = np.zeros((350,l2))
	Pmod_dt_prime = np.zeros((350,l2))
	Pmod_tt_prime = np.zeros((350,l2))
	
	for count,iz in enumerate(red2):
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/A_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Aprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/B_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Bprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/C_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Cprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/D_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Dprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/E_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Eprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/F_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Fprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/G_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Gprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/H_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Hprime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/Pmod_dd_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Pmod_dd_prime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/Pmod_dt_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Pmod_dt_prime[:,count] = f[:,1]
		#------------------------
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
		'/PT_coeff/Pmod_tt_'+str(iz)+'.txt')
		f = np.loadtxt(dat_file_path)
		kpt = f[:,0]
		Pmod_tt_prime[:,count] = f[:,1]
			
	### interpolate the pt coeff on the chosen scale
	A = np.zeros((len(kclass),l2))
	B = np.zeros((len(kclass),l2))
	C = np.zeros((len(kclass),l2))
	D = np.zeros((len(kclass),l2))
	E = np.zeros((len(kclass),l2))
	F = np.zeros((len(kclass),l2))
	G = np.zeros((len(kclass),l2))
	H = np.zeros((len(kclass),l2))
	Pmod_dd = np.zeros((len(kclass),l2))
	Pmod_dt = np.zeros((len(kclass),l2))
	Pmod_tt = np.zeros((len(kclass),l2))
	for i in xrange(l2):
		A[:,i] = np.interp(kclass, kpt, Aprime[:,i]) 
		B[:,i] = np.interp(kclass, kpt, Bprime[:,i]) 
		C[:,i] = np.interp(kclass, kpt, Cprime[:,i]) 
		D[:,i] = np.interp(kclass, kpt, Dprime[:,i]) 
		E[:,i] = np.interp(kclass, kpt, Eprime[:,i]) 
		F[:,i] = np.interp(kclass, kpt, Fprime[:,i]) 
		G[:,i] = np.interp(kclass, kpt, Gprime[:,i]) 
		H[:,i] = np.interp(kclass, kpt, Hprime[:,i]) 
		Pmod_dd[:,i] = np.interp(kclass, kpt, Pmod_dd_prime[:,i]) 
		Pmod_dt[:,i] = np.interp(kclass, kpt, Pmod_dt_prime[:,i]) 
		Pmod_tt[:,i] = np.interp(kclass, kpt, Pmod_tt_prime[:,i]) 
		
		
		
	return A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt
