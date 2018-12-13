
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
import error
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys


def tuning(self, kclass, Phh, redshift, znumber, case, Massbins, dim):
	####################################################################
	#### select the k mode according to the ones in Raccanelli et al. 2017
	k_article = [0.35,0.2,0.15]
	
	if case == 1:
		kmax = k_article[0]
	elif case == 2:
		kmax = k_article[1]
	elif case == 3:
		kmax = k_article[2]
		
	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	####################################################################
	#### create a scale array limited by kmin and kmax
	try:
		self.kmax
	except:
		self.kmax = False
	if self.kmax and self.kmax < kmax:
		lim_h = np.where(kclass <= self.kmax)[0]
	else:
		lim_h = np.where(kclass <= kmax)[0]

	klim_h = np.amax(lim_h) # define the higher limit
	
	#-------------------------------------------------------------------
	Vs = 1000**3 # for a periodic box of 1000 h-1Mpc
	kmin = 2 * math.pi * Vs**(-1/3.)
	
	try:
		self.kmin
	except:
		self.kmin = False

	if self.kmin and self.kmin > kmin:
		lim_l = np.where(kclass >= self.kmin)[0]
	else:
		lim_l = np.where(kclass >= kmin)[0]
	#------------------------------------------------------------------
	##### =====> 
	kclass = kclass[lim_l[0]:klim_h+1]
	Phh = Phh[lim_l[0]:klim_h+1]

	####################################################################
	#### interpolate on selected redshift
	
	if dim == '2d':
		Phhbis = np.zeros((len(kclass),znumber))
		for ik in xrange(len(kclass)):
			f = interp1d(red2, Phh[ik,:], kind='cubic', fill_value='extrapolate')
			Phhbis[ik,:] = f(redshift)
	if dim == '3d':
		Phhbis = np.zeros((len(kclass),znumber,len(Massbins)))
		for j in xrange(len(Massbins)):
			for ik in xrange(len(kclass)):
				f = interp1d(red2, Phh[ik,:,j], kind='cubic', fill_value='extrapolate')
				Phhbis[ik,:,j] = f(redshift)
			
	return kclass, Phhbis
