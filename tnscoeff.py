

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




def tnscoeff(self, data, case, kclass, Massbins):
       
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

	AB2 = np.zeros((350,l2,len(Massbins)))
	AB4 = np.zeros((350,l2,len(Massbins)))
	AB6 = np.zeros((350,l2,len(Massbins)))
	AB8 = np.zeros((350,l2,len(Massbins)))
	
	if case == 'lin':
		for count,iz in enumerate(red2):
			for count2,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				if ind2 == 0:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_lin1_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 1:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_lin2_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 2:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_lin3_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 3:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_lin4_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
	#-------------------------------------------------------------------
	if case == 'pl':
		for count,iz in enumerate(red2):
			for count2,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				if ind2 == 0:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pl1_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 1:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pl2_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 2:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pl3_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 3:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pl4_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
	#-------------------------------------------------------------------
	if case == 'exp':
		for count,iz in enumerate(red2):
			for count2,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				if ind2 == 0:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pt1_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 1:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pt2_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 2:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pt3_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
				elif ind2 == 3:
					dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
					'/TNS_coeff/tns_pt4_z='+str(iz)+'.txt')
					f = np.loadtxt(dat_file_path)
					AB2[:,iz,count2] = f[:,0]
					AB4[:,iz,count2] = f[:,1]
					AB6[:,iz,count2] = f[:,2]
					AB8[:,iz,count2] = f[:,3]
		
			
	### interpolate the pt coeff on the chosen scale
	
	dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV'\
	'/PT_coeff/A_0.0.txt')
	f = np.loadtxt(dat_file_path)
	kpt = f[:,0]
	AB2bis = np.zeros((len(kclass),l2,len(Massbins)))
	AB4bis = np.zeros((len(kclass),l2,len(Massbins)))
	AB6bis = np.zeros((len(kclass),l2,len(Massbins)))
	AB8bis = np.zeros((len(kclass),l2,len(Massbins)))

	for m in xrange(len(Massbins)):
		for i in xrange(l2):
			AB2bis[:,i,m] = np.interp(kclass, kpt, AB2[:,i,m]) 
			AB4bis[:,i,m] = np.interp(kclass, kpt, AB4[:,i,m]) 
			AB6bis[:,i,m] = np.interp(kclass, kpt, AB6[:,i,m]) 
			AB8bis[:,i,m] = np.interp(kclass, kpt, AB8[:,i,m]) 

	return AB2bis, AB4bis, AB6bis, AB8bis
