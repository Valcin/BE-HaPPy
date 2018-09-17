
from classy import Class
from matplotlib.colors import LogNorm
from bcoeff import bcoeff
from ls_coeff import lscoeff
from pt_coeff import ptcoeff
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import os
import numpy as np
import warnings
import csv
import sys
sys.path.append('/home/david/codes/FAST-PT')
import myFASTPT as FPT
from bias import Halo


def ci(self, cosmo, data, red2):

####################################################################
	#### import the requested redshift(s) 
	try:
		self.z
	except:
		self.z = []

	if len(self.z)>0:
		redshift = self.z
	#--------------------------------------------
	try:
		self.redshift
	except:
		self.redshift = False

	if self.redshift:
		redshift = self.redshift
	#--------------------------------------------
	if not len(self.z)>0 and not self.redshift:
		raise ValueError('Please define redshift(s) named redshift or z')
	
	####################################################################
	#### Check if the total neutrino mass corresponds to one of the available ones
	#### get_current_derived_parameters returns a dict so must be converted
	m = cosmo.get_current_derived_parameters(['m_ncdm_tot'])
	m = m.values()
	m = [ round(elem, 2) for elem in m ]
	mv = [0.0, 0.03, 0.06, 0.10, 0.13, 0.15, 0.30]
	if m[0] not in mv:
		raise ValueError('Sorry the code is only available for Mv = 0.0, 0.03, 0.06, 0.10, 0.13, 0.15, 0.30 and your Mv is '+str(m[0])+'. Please modify you total neutrino mass.')

	

	####################################################################
	#### Store the selected redshifts in a array and deduce its length for the loops
	#### array manipulation because len() and size only work for znumber >1
	a = np.array(redshift)
	znumber = a.size 
	redshift = np.zeros(znumber,'float64') 
	redshift[:] = a

	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	#~ red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	
	####################################################################
	#### Store the mass bins available in an array
	mbins = ['M1', 'M2', 'M3', 'M4']

	####################################################################
	#### Since the cosmo.pk k's are bounded in [0.000000e+00:5.366287e+00]
	#### we must extract the k values from the get_transfer list so they coincide. 
	kget = cosmo.get_transfer(red2[0])
	kclass = kget.get('k (h/Mpc)')
	bound = np.where((kclass > 0)&(kclass < 5.366287e+00))[0]
	kclass = kclass[bound]

		
	####################################################################
	#### get the value of h to rescale the power spectrum and wavenumber
	#### because of the difference between Class and Class python wrapper 
	param = cosmo.get_current_derived_parameters(['h','Omega0_lambda'])
	param = param.values()
	param = [ round(elem, 5) for elem in param ]
	h = param[0]
	Omega_lambda = param[1]

	####################################################################
	#### get the transfer function from class
	d_b = np.zeros((len(kclass), l2), 'float64')
	d_cdm = np.zeros((len(kclass), l2), 'float64')
	d_tot = np.zeros((len(kclass), l2), 'float64')
	for i in xrange(l2):
		transfer = cosmo.get_transfer(red2[i])
		d_b[:,i] = transfer.get('d_b')[bound]
		d_cdm[:,i] = transfer.get('d_cdm')[bound]
		d_tot[:,i] = transfer.get('d_tot')[bound]

	####################################################################
	#### import Omega_b and Omega_cdm from class. Remember to add Omega_cdm in classy and recompile after
	Omega_b = cosmo.Omega_b()
	Omega_cdm = cosmo.Omega_cdm()
	Omega_m = cosmo.Omega_m()


	####################################################################
	#### define the CDM + baryons transfer function 
	T_cb = np.zeros((len(kclass), l2), 'float64')
	T_cb = (Omega_cdm * d_cdm + Omega_b * d_b)/(Omega_cdm + Omega_b)
    
    
	####################################################################
	#### get the non linear power spectrum from class
	pk = np.zeros((len(kclass), l2), 'float64')
	for ik in xrange(len(kclass)):
		for iz in xrange(l2):
			pk[ik,iz] = cosmo.pk(kclass[ik], red2[iz])

	
	####################################################################
	#### get the linear power spectrum from class
	#~ pk_lin = np.zeros((len(kclass), l2), 'float64')
	#~ for ik in xrange(len(kclass)):
		#~ for iz in xrange(l2):
			#~ pk_lin[ik,iz] = cosmo.pk_lin(kclass[ik], red2[iz])	
			
	####################################################################
	#### Define the linear growth factor and growth rate (growth factor f in class)
	f = np.zeros(znumber, 'float64')
	D = np.zeros(znumber, 'float64')
	for iz in xrange(znumber):
		f[iz] = cosmo.scale_independent_growth_factor_f(redshift[iz])
		D[iz] = cosmo.scale_independent_growth_factor(redshift[iz])
		
	
	return redshift, m[0], h, d_tot, T_cb, pk, f, D
