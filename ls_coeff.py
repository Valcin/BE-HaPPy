
from classy import Class
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from bcoeff import bcoeff
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys


def lscoeff(self, data, mv, Massbins):

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

	####################################################################
	#### get the rescaling coefficients according to neutrino mass
	#### if you want to change the large sclae bias from Tinker to Crocce + ST 
	#### change the line 'LS_z='+str(i)+'_.txt' in 'LS2_z='+str(i)+'_.txt'
	bcc_LS000 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.0eV/large_scale/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS000[ind,count] = f[ind2]
	#------------------------------
	bcc_LS003 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.03/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS003[ind,count] = f[ind2]
	#------------------------------
	bcc_LS006 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.06/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS006[ind,count] = f[ind2]
	#------------------------------
	bcc_LS010 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.10/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS010[ind,count] = f[ind2]
	#------------------------------
	bcc_LS013 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.13/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS013[ind,count] = f[ind2]
	#------------------------------
	bcc_LS015 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/0.15eV/large_scale/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS015[ind,count] = f[ind2]
	#------------------------------
	bcc_LS030 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.30/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS030[ind,count] = f[ind2]

	###############################################################
		if mv == 0.0:
			return bcc_LS000, bcc_LS000
		if mv == 0.03:
			return bcc_LS000, bcc_LS003
		elif mv == 0.06:
			return bcc_LS000, bcc_LS006
		elif mv == 0.10:
			return bcc_LS000, bcc_LS010
		elif mv == 0.13:
			return bcc_LS000, bcc_LS013
		elif mv == 0.15:
			return bcc_LS000, bcc_LS015
		elif mv == 0.30:
			return bcc_LS000, bcc_LS030
		
		
