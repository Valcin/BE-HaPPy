
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
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.1/'\
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
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.3/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS030[ind,count] = f[ind2]
	#------------------------------
	bcc_LS045 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.45/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS045[ind,count] = f[ind2]
	#------------------------------
	bcc_LS060 = np.zeros((l2,len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'BE_HaPPy/coefficients/other neutrinos masses/0.6/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS060[ind,count] = f[ind2]

	###############################################################
	av_mass = [0.0,0.03,0.06,0.10,0.13,0.15,0.30,0.45,0.60]
	dim = np.shape(bcc_LS000)
	bcc_LSmassive = np.zeros((dim[0],dim[1]))
	for m in xrange(0, dim[0]):
		for n in xrange(0, dim[1]):
			bvalue = [bcc_LS000[m,n],bcc_LS003[m,n],bcc_LS006[m,n],bcc_LS010[m,n],bcc_LS013[m,n],bcc_LS015[m,n],\
			bcc_LS030[m,n],bcc_LS045[m,n],bcc_LS060[m,n]]
			fLS = interp1d(av_mass, bvalue)
			bcc_LSmassive[m,n] = fLS(mv)
			
	return bcc_LS000, bcc_LSmassive


