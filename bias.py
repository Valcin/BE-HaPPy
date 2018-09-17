
from classy import Class
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from cimp import cimp
from bcoeff import bcoeff
from ls_coeff import lscoeff
from pt_coeff import ptcoeff
import error
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys


############################## INPUT ###################################
# neutrino parameters
#hierarchy = 'degenerate' #'degenerate', 'normal', 'inverted'
#Mnu       = 0.0, 0.03, 0.06, 0.10, 0.13, 0.15, 0.30  #eV
#Nnu       = 0, 3 #number of massive neutrinos
#Neff      = 3.046

# cosmological parameters
#h       = 0.6711
#Omega_c = 0.2685 - Mnu/(93.14*h**2)
#Omega_b = 0.049
#Omega_k = 0.0
#tau     = None

# initial P(k) parameters
#ns           = 0.9624
#As           = 2.13e-9
#pivot_scalar = 0.05
#pivot_tensor = 0.05
########################################################################

def Halo(self, cosmo, data, model, case, Massbins, err = None):

       
	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	####################################################################
	#### Store the mass bins available in an array
	mbins = ['M1', 'M2', 'M3', 'M4']

	
	####################################################################
	#### get the coefficients from the dat files in the data directory 

	#-------------------------------------------------------------------
	if model =='exp':
		b1, b2, bs, b3nl = bcoeff(self, data, model, case, Massbins)

	#-------------------------------------------------------------------
	if model =='pl':
		b1, b2, b3, b4 = bcoeff(self, data, model, case, Massbins)
	
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
	#### import classy results
	redshift, znumber, mv, h, d_tot, T_cb, kclass, pk, f, D = cimp(self, cosmo, data)
	
	####################################################################
	###### compute the power spectrum with linear bias


	if model == 'lin':
		# tinker effective bias
		bcc = lscoeff(self,data, mv,Massbins)[1]
					
		# compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.
		bmm = np.zeros((len(kclass),l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				bmm[:,iz, count] = bcc[iz,count] * (T_cb[:,iz]/d_tot[:,iz])# * (bcc_LS0[iz]/denom[iz])

		# Compute the halo Power spectrum in real space
		Phh = np.zeros((len(kclass), l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				Phh[:,iz,count] = pk[:,iz] * bmm[:,iz, count]**2
			
			
		# rescale the k and power spectrum because of classy/class difference
		kclass /= h
		Phh *= h**3
		
		# create a scale array limited by kmin and kmax
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
		
		# interpolate on selected redshift
		Phhbis = np.zeros((len(kclass),znumber,len(Massbins)))
		for j in xrange(len(Massbins)):
			for ik in xrange(len(kclass)):
				f = interp1d(red2, Phh[ik,:,j], kind='cubic', fill_value='extrapolate')
				Phhbis[ik,:,j] = f(redshift)
				
		if err == True:
			if mv == 0.0 or mv == 0.15:
				error.error(self, data, kclass, Phhbis, redshift, mv, Massbins)
			else:
				raise ValueError('the simulation spectra are only available for Mv = 0.0 or 0.15eV sorry')
		
		return kclass, Phhbis
		
	####################################################################
	###### compute the bias and halo power spectrum for power law model


	elif model == 'pl':
		bcc = np.zeros((len(kclass), l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bcc[:,iz, count] = b1[iz,count] + b2[iz,count]*(kclass**2) + b3[iz,count]*(kclass**3) \
				+ b4[iz,count]*(kclass**4) 
				
		bcc_LS000 = lscoeff(self,data, mv,Massbins)[0]
		bcc_LSmassive = lscoeff(self,data, mv,Massbins)[1]
		# if mv = 0.0eV bcc_LS000 = bcc_LSmassive
		for iz in xrange(l2):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LSmassive[iz,count]/bcc_LS000[iz,count]
		
				
		# compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.
		bmm = np.zeros((len(kclass), l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bmm[:,iz, count] = bcc[:,iz,count] * (T_cb[:,iz]/d_tot[:,iz])

		# Compute the halo Power spectrum in real space
		Phh = np.zeros((len(kclass), l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				Phh[:,iz,count] = pk[:,iz] * bmm[:,iz, count]**2
			
			
		# rescale the k and power spectrum because of classy/class difference
		kclass /= h
		Phh *= h**3
		
		# create a scale array limited by kmin and kmax
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
		
		#~ # interpolate on selected redshift
		Phhbis = np.zeros((len(kclass),znumber,len(Massbins)))
		for j in xrange(len(Massbins)):
			for ik in xrange(len(kclass)):
				f = interp1d(red2, Phh[ik,:,j], kind='cubic', fill_value='extrapolate')
				Phhbis[ik,:,j] = f(redshift)
				
		if err == True:
			if mv == 0.0 or mv == 0.15:
				error.error(self, data, kclass, Phhbis, redshift, mv, Massbins)
			else:
				raise ValueError('the simulation spectra are only available for Mv = 0.0 or 0.15eV sorry')
		
		return kclass, Phhbis
		
		
	####################################################################
	###### compute the one loop correction with FAST-PT for the expansion model


	elif model == 'exp':
		
		A,B,C,D,E,F,G,H,Pmod_dd,_,_ = ptcoeff(self, data, kclass, Massbins)		
		
		# compute the halo power spectrum given the coefficient
		PhhDD = np.zeros((len(kclass),l2,len(Massbins)))
		PhhDT = np.zeros((len(kclass),l2,len(Massbins)))
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				# density spectrum rescaled by transfer function from bcc ----> bmm
				PhhDD[:,iz,count] = b1[iz,count]**2*Pmod_dd[:,iz] + b1[iz,count]*b2[iz,count]*A[:,iz] + 1/4.*b2[iz,count]**2*B[:,iz] + \
				b1[iz,count]*bs[iz,count]*C[:,iz] + 1/2.*b2[iz,count]*bs[iz,count]*D[:,iz] + 1/4.*bs[iz,count]**2*E[:,iz] +\
				2*b1[iz,count]*b3nl[iz,count]*F[:,iz] #* (T_cb[:,iz]/d_tot[:,iz])**2
				# cross velocity spectrum
				#~ PhhDT[:,iz,count] = b1[iz,count]* Pmod_dt[:,iz] + b2[iz,count]*G[:,iz] + bs[iz,count]*H[:,iz] + b3nl[iz,count]*F[:,iz] \
				#~ *(T_cb[:,iz]/d_tot[:,iz])
		
		bcc_LS000 = lscoeff(self,data, mv,Massbins)[0]
		bcc_LSmassive = lscoeff(self,data, mv,Massbins)[1]
		# if mv = 0.0eV bcc_LS000 = bcc_LSmassive
		for iz in xrange(l2):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LSmassive[iz,count]/bcc_LS000[iz,count]
		
		# create a scale array limited by kmin and kmax

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
		PhhDD = PhhDD[lim_l[0]:klim_h+1]
		
		# interpolate on selected redshift
		PhhDDbis = np.zeros((len(kclass),znumber,len(Massbins)))
		for j in xrange(len(Massbins)):
			for ik in xrange(len(kclass)):
				f = interp1d(red2, PhhDD[ik,:,j], kind='cubic', fill_value='extrapolate')
				PhhDDbis[ik,:,j] = f(redshift)
				
		if err == True:
			if mv == 0.0 or mv == 0.15:
				error.error(self, data, kclass, PhhDDbis, redshift, mv, Massbins)
			else:
				raise ValueError('the simulation spectra are only available for Mv = 0.0 or 0.15eV sorry')

		return kclass, PhhDDbis
		
		
	





