from montepython.likelihood_class import Likelihood
import os
import numpy as np
import warnings
from numpy import newaxis as na
from math import exp, log, pi, log10
#~ sys.path.append('/home/david/codes/BE_HaPPy')
#~ from BE_HaPPy.bias import Halo
#~ from BE_HaPPy.red_ps import rspec
import time
import math
import sys
sys.path.append('/home/david/codes/FAST-PT')
import myFASTPT as FPT

class BE_HaPPy(Likelihood):
	def __init__(self, path, data, command_line):

		Likelihood.__init__(self, path, data, command_line)

		self.need_cosmo_arguments(data, {'output': 'mPk mTk'})
		self.need_cosmo_arguments(data, {'z_max_pk': self.zmax}) # needed for transfer
		self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': 1.9*self.kmax})
		self.need_cosmo_arguments(data, {'non linear': 'halofit'})


		#~ THINK OF HST PRIORS
		# Else the file will be created in the loglkl() function.
		return

	def loglkl(self, cosmo, data):
		start = time.time()
		####################################################################
		####################################################################
		# Convenience variables: store the nuisance parameters in short named
		# variables
		A_shot = (data.mcmc_parameters['A_shot']['current'] *
		data.mcmc_parameters['A_shot']['scale'])
		sigma_v = (data.mcmc_parameters['sigma_v']['current'] *
		data.mcmc_parameters['sigma_v']['scale'])
		b1 = (data.mcmc_parameters['b1']['current'] *
		data.mcmc_parameters['b1']['scale'])
		b2 = (data.mcmc_parameters['b2']['current'] *
		data.mcmc_parameters['b2']['scale'])
		
		# compute bs and b3nl according to local lagrangian bias
		bs = -4/7. * (b1 -1)
		b3nl = 32/315. * (b1 -1)
		
		####################################################################
		####################################################################
		#### import classy results
		#### import the requested redshift(s) 
		try:
			self.z
		except:
			self.z = False

		#~ if len(self.z)>0:
		if self.z:
			redshift = self.z
		#--------------------------------------------
		try:
			self.redshift
		except:
			self.redshift = False

		if self.redshift:
			redshift = self.redshift
		#--------------------------------------------
		if not self.z and not self.redshift:
			raise ValueError('Please define redshift(s) named redshift or z')
		
		#### Check if the total neutrino mass corresponds to one of the available ones
		#### get_current_derived_parameters returns a dict so must be converted
		m = cosmo.get_current_derived_parameters(['m_ncdm_tot'])
		m = m.values()
		m = [ round(elem, 2) for elem in m ]
		if m[0] < 0 or m[0] > 0.60:
			raise ValueError('Sorry the code is only available for Mv = [0.0 - 0.60] your Mv is '+str(m[0])+'. Please modify you total neutrino mass.')

		#~ print 'The total neutrino mass is '+str(m[0])

		#### Store the selected redshifts in a array and deduce its length for the loops
		#### array manipulation because len() and size only work for znumber >1
		a = np.array(redshift)
		znumber = a.size 
		redshift = np.zeros(znumber,'float64') 
		redshift[:] = a

		#### Store the redshifts where bcc fit  and bcc Ls are available in arrays
		red2 = [0.0,0.5,1.0,2.0]
		l2= len(red2)
		
		#### get the value of h to rescale the power spectrum and wavenumber
		#### because of the difference between Class and Class python wrapper 
		param = cosmo.get_current_derived_parameters(['h','Omega0_lambda'])
		param = param.values()
		param = [ round(elem, 5) for elem in param ]
		h = param[0]
		Omega_lambda = param[1]

		#### get the transfer function from class
		kget = cosmo.get_transfer(self.z)
		kclass = kget.get('k (h/Mpc)')
		d_b = np.zeros((len(kclass)), 'float64')
		d_cdm = np.zeros((len(kclass)), 'float64')
		d_tot = np.zeros((len(kclass)), 'float64')
		transfer = cosmo.get_transfer(self.z)
		d_b = transfer.get('d_b')
		d_cdm = transfer.get('d_cdm')
		d_tot = transfer.get('d_tot')
			
		#### Since the cosmo.pk k's are bounded in [0.000000e+00:5.366287e+00]
		#### and Fast- PT requires a evenly sampled array in log scale 
		#### we interpolate
		kbound = np.logspace(np.log10(self.kmin), np.log10(self.kmax), self.kbins)
		d_b = np.interp(kbound, kclass,d_b)
		d_cdm = np.interp(kbound, kclass,d_cdm)
		d_tot = np.interp(kbound, kclass,d_tot)

		#### import Omega_b and Omega_cdm from class. Remember to add Omega_cdm in classy and recompile after
		Omega_b = cosmo.Omega_b()
		Omega_cdm = cosmo.Omega_cdm()
		Omega_m = cosmo.Omega_m()

		#### define the CDM + baryons transfer function 
		T_cb = np.zeros((len(kclass)), 'float64')
		T_cb = (Omega_cdm * d_cdm + Omega_b * d_b)/(Omega_cdm + Omega_b)
		
		#### get the linear power spectrum from class
		pk_lin = np.zeros((len(kbound)), 'float64')
		for ik in xrange(len(kbound)):
			pk_lin[ik] = cosmo.pk_lin(kbound[ik], self.z)
				
		#### get the non linear power spectrum from class
		pk = np.zeros((len(kbound)), 'float64')
		for ik in xrange(len(kbound)):
			pk[ik] = cosmo.pk(kbound[ik], self.z)
				
		#### Define the linear growth factor and growth rate (growth factor f in class)
		f = Omega_m**0.55
		D = cosmo.scale_independent_growth_factor(self.z)
		
		####################################################################
		####################################################################
		### compute perturbation terms 
		
		# set the parameters for the power spectrum window and
		# Fourier coefficient window 
		#P_window=np.array([.2,.2])  
		C_window=0.95

		# padding length 
		nu=-2; n_pad=len(kbound)
		n_pad=int(0.5*len(kbound))
		to_do=['all']
						
		# initialize the FASTPT class 
		# including extrapolation to higher and lower k  
		# time the operation
		fastpt=FPT.FASTPT(kbound,to_do=to_do,n_pad=n_pad, verbose=True) 
			
		# calculate 1loop SPT (and time the operation) for density
		P_spt_dd=fastpt.one_loop_dd(pk_lin,C_window=C_window)
			
		# calculate 1loop SPT (and time the operation) for velocity
		P_spt_tt=fastpt.one_loop_tt(pk_lin,C_window=C_window)
			
		# calculate 1loop SPT (and time the operation) for velocity - density
		P_spt_dt=fastpt.one_loop_dt(pk_lin,C_window=C_window)
			
			
		#calculate tidal torque EE and BB P(k)
		#~ P_RSD=fastpt.RSD_components(P,1.0,C_window=C_window)	

		# update the power spectrum
		Pmod_dd=pk_lin+P_spt_dd[0]
		Pmod_dt=pk_lin+P_spt_dt[0]
		Pmod_tt=pk_lin+P_spt_tt[0]	
		A = P_spt_dd[2]
		B = P_spt_dd[3]
		C = P_spt_dd[4]
		D = P_spt_dd[5]
		E = P_spt_dd[6]
		F = P_spt_dd[7]
		G = P_spt_dt[2]
		H = P_spt_dt[3]

		####################################################################
		####################################################################
		# compute the halo power spectrum given the coefficient
		PhhDD = np.zeros((len(kbound)))
		# density spectrum rescaled by transfer function from bcc ----> bmm
		PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*bs*C + 1/2.*b2*bs*D + 1/4.*bs**2*E +\
		2*b1*b3nl*F 
	
		# cross velocity spectrum
		PhhDT = np.zeros((len(kbound)))
		PhhDT = b1* Pmod_dt + b2*G + bs*H + b3nl*F 
			
		
		# compute tns coeff  and interpolate on z	
		AB2, AB4, AB6, AB8 = fastpt.RSD_ABsum_components(pk_lin,f,b1,C_window=C_window)

		# get the selected k array
		kappa = np.zeros((len(kbound)))
		coeffA = np.zeros((len(kbound)))
		coeffB = np.zeros((len(kbound)))
		coeffC = np.zeros((len(kbound)))
		coeffD = np.zeros((len(kbound)))
		coeffE = np.zeros((len(kbound)))

		kappa = kbound*(data.mcmc_parameters['sigma_v']['current']*data.mcmc_parameters['sigma_v']['scale'])*f*D
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		
		# compute ps in redshift space
		Pred = np.zeros((len(kbound)))
		Pred = PhhDD*coeffA  + 2/3.*f*PhhDT*coeffB + 1/5.*f**2*Pmod_tt*coeffC + 1/3.*AB2*coeffB \
		+ 1/5.*AB4*coeffC + 1/7.*AB6*coeffD + 1/9.*AB8*coeffE

		end = time.time()
		print 'total time is '+str((end - start))
		return -0.5 * chi2
