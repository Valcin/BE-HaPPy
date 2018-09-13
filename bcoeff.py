
from classy import Class
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys




def bcoeff(self, data, model, case, Massbins):
	
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

	#-------------------------------------------------------------------
	if model =='exp':
	
		b1 = np.zeros((l2,len(Massbins)))
		b2 = np.zeros((l2,len(Massbins)))
		bs = np.zeros((l2,len(Massbins)))
		b3nl = np.zeros((l2,len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case1/coeff_3exp_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					bs[ind,count] = f[ind2,2]
					b3nl[ind,count] = f[ind2,3]

		if case == 2:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_3exp_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					bs[ind,count] = f[ind2,2]
					b3nl[ind,count] = f[ind2,3]

		if case == 3:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_3exp_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					bs[ind,count] = f[ind2,2]
					b3nl[ind,count] = f[ind2,3]
					
		return b1, b2, bs, b3nl

	#-------------------------------------------------------------------
	if model =='pl':
	
		b1 = np.zeros((l2,len(Massbins)))
		b2 = np.zeros((l2,len(Massbins)))
		b3 = np.zeros((l2,len(Massbins)))
		b4 = np.zeros((l2,len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case1/coeff_pl_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					b3[ind,count] = f[ind2,2]
					b4[ind,count] = f[ind2,3]

		if case == 2:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_pl_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					b3[ind,count] = f[ind2,2]
					b4[ind,count] = f[ind2,3]

		if case == 3:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_pl_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					b3[ind,count] = f[ind2,2]
					b4[ind,count] = f[ind2,3]

		return b1, b2, b3, b4

	
