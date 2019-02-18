from montepython.likelihood_class import Likelihood
import os
import numpy as np
import warnings
from numpy import newaxis as na
from math import exp, log, pi, log10
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

		#~ #-------------------------------------------------
		#~ #---------------- Data ---------------------------
		#~ #-------------------------------------------------
		Class = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/class/test_z2_pk.dat')
		self.kclass = Class[:,0]
		self.Pclass = Class[:,1]
		
		
		# store the calibrated redshift
		self.red = [0.0,0.5,1.0,2.0] 
		self.lred = len(self.red) 
		ind = self.red.index(self.z)
		
		print 'Total neutrino mass is ' + str(self.Mnu)+'eV'
		print 'you chose the mass bin M' + str(self.mbin +1)
		if self.rsd == 1:
				print 'you chose the Kaiser model'
		elif self.rsd == 2:
				print 'you chose the Scoccimaro model'
		elif self.rsd == 3:
				print 'you chose the TNS model'
		if self.bmodel == 1:
			print 'you chose the linear bias'
		elif self.bmodel == 2:
			print 'you chose the polynomial bias'
		elif self.bmodel == 3:
			print 'you chose the perturbation theory bias'
		print ('')
		
		#### load halo redshift space ps
		### for test
		#~ if err == True:
			#~ if mv == 0.0 or mv == 0.15:
				#~ error.error(self, data, kclass, Phhbis, redshift, mv, Massbins)
			#~ else:
				#~ raise ValueError('the simulation spectra are only available for Mv = 0.0 or 0.15eV sorry')
		#~ if self.Mnu != 0.0 or self.Mnu != 0.15:
				#~ raise ValueError('the simulation spectra are only available for Mv = 0.0 or 0.15eV sorry')
				self.data_directory = data.path['root']
		print self.data_directory
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/simu_ps')
		sys.path.append(dat_file_path)
		from load_data import ld_data
		kcamb, Pcamb, k, Pmm, PH1, PH2, PH3 , PH4, errPhh1, errPhh2, errPhh3, errPhh4, bias1, bias2, bias3, bias4, \
		bias1s, bias2s, bias3s, bias4s, errb1, errb2, errb3, errb4, Pmono1, Pmono2, Pmono3, Pmono4, errPr1, errPr2, errPr3,\
		errPr4, kclass, Tm, Tcb, noise1, noise2, noise3, noise4 = ld_data(self.Mnu, self.red, ind)
		
		self.ksimu = k
		if self.mbin == 0:
			self.Psimu = Pmono1
			self.err = errPr1
		elif self.mbin == 1:
			self.Psimu = Pmono2
			self.err = errPr2
		elif self.mbin == 2:
			self.Psimu = Pmono3
			self.err = errPr3
		elif self.mbin == 3:
			self.Psimu = Pmono4
			self.err = errPr4
		
		return

	def loglkl(self, cosmo, data):
		start = time.time()
		####################################################################
		####################################################################
		# Convenience variables: store the nuisance parameters in short named
		# variables
		A_shot = (data.mcmc_parameters['A_shot']['current'] *
		data.mcmc_parameters['A_shot']['scale'])
		
		
		####################################################################
		####################################################################
		#### import classy results
		#### import the requested redshift(s) 
		if self.z == 0:
			redshift = self.z
		else:
			try:
				self.z
			except:
				self.z = False

			if self.z:
				redshift = self.z
				
			if not self.z:
				raise ValueError('Please define redshift(s)')

		#### Store the selected redshifts in a array and deduce its length for the loops
		#### array manipulation because len() and size only work for znumber >1
		red = np.array(redshift)
		znumber = red.size 
		redshift = np.zeros(znumber,'float64') 
		redshift[:] = red

		#### get the transfer function from class
		#~ kget = cosmo.get_transfer(self.z, output_format='camb')
		#~ kclass = kget.get('k (h/Mpc)')
		#~ d_b = np.zeros((len(kclass)), 'float64')
		#~ d_cdm = np.zeros((len(kclass)), 'float64')
		#~ d_tot = np.zeros((len(kclass)), 'float64')
		#~ transfer = cosmo.get_transfer(self.z)
		#~ d_b = transfer.get('d_b')
		#~ d_cdm = transfer.get('d_cdm')
		#~ d_tot = transfer.get('d_tot')
		
		#### import Omega_b and Omega_cdm from class. Remember to add Omega_cdm in classy and recompile after
		Omega_b = cosmo.Omega_b()
		Omega_cdm = cosmo.Omega_cdm()
		Omega_m = cosmo.Omega_m()
		Omega_Lambda = cosmo.Omega_Lambda()
		Omega_k = cosmo.Omega0_k()
		h = cosmo.h()
		
		#~ print h, Omega_m, Omega_b, Omega_cdm, Omega_Lambda, Omega_k
		
		# compute the scale array from the data file
		kbound = np.logspace(np.log10(self.kmin), np.log10(self.kmax), self.kbins)
		
		
		#### get the linear power spectrum from class. here multiply input k array by h because get_pk uses 1/mpc 
		if self.cdm == 0:
			pk_lin = cosmo.get_pk_array(kbound*h, redshift, len(kbound), znumber, 0) #if we want Pmm
		elif self.cdm == 1:
			pk_lin = cosmo.get_pk_cb_array(kbound*h, redshift, len(kbound), znumber, 0) # if we want Pcb
			
		
		### compare classy amplitude with classy
		import matplotlib.pyplot as plt
		#~ plt.plot(self.kclass, self.Pclass, c='r')
		#~ plt.plot(kbound, pk_lin*h**3, c='b')
		#~ plt.xscale('log')
		#~ plt.yscale('log')
		#~ plt.show()
		
		### rescale the amplitude of pk_lin accoridngly
		pk_lin *= h**3 
				
		#### Define the linear growth factor and growth rate (growth factor f in class)
		#~ fz = Omega_m**0.55
		fz = cosmo.scale_independent_growth_factor_f(self.z)
		Dz = cosmo.scale_independent_growth_factor(self.z)
		#~ print fz, Dz
		
		####################################################################
		####################################################################
		### load perturbation terms 
		pt_terms = ['Pmod_dd', 'Pmod_dt', 'Pmod_tt','A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
		
		
		Pmod_dd = self.interpol_pt(pt_terms[0], kbound, self.z)
		Pmod_dt = self.interpol_pt(pt_terms[1], kbound, self.z)
		Pmod_tt = self.interpol_pt(pt_terms[2], kbound, self.z)
		A = self.interpol_pt(pt_terms[3], kbound, self.z)
		B = self.interpol_pt(pt_terms[4], kbound, self.z)
		C = self.interpol_pt(pt_terms[5], kbound, self.z)
		D = self.interpol_pt(pt_terms[6], kbound, self.z)
		E = self.interpol_pt(pt_terms[7], kbound, self.z)
		F = self.interpol_pt(pt_terms[8], kbound, self.z)
		G = self.interpol_pt(pt_terms[9], kbound, self.z)
		H = self.interpol_pt(pt_terms[10], kbound, self.z)
		
		
		####################################################################
		####################################################################
		### load bias coefficients
		
		b1, b2, b3, b4 = self.bcoeff(self.z)
		
		####################################################################
		####################################################################
		### load bias coefficients

		alpha  = self.rescaling()
		
		
		####################################################################
		####################################################################
		### compute the redshift power spectrum
		
		Pred = self.red_ps(data, kbound, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt, alpha)

		####################################################################
		####################################################################

		### interpolate data on kbound
		if len(self.Psimu) != len(kbound):
			self.Psimu = np.interp(kbound, self.ksimu, self.Psimu)
			self.err = np.interp(kbound, self.ksimu, self.err)
		
		### plot to test
		#~ plt.plot(kbound, self.Psimu, c='b')
		#~ plt.plot(kbound, Pred, c='r')
		#~ plt.xscale('log')
		#~ plt.yscale('log')
		#~ plt.show()
		
		
		### compute the chi square
		inv_sigma2 = 1.0/(self.err**2)
		chi2 = -0.5*(np.sum((self.Psimu-Pred)**2*inv_sigma2 - np.log(inv_sigma2)))
			
		end = time.time()
		print 'total time is '+str((end - start))
		return -chi2
		
		
##########################################################################################################################
##########################################################################################################################
		
	def interpol_pt(self, pt_term, kbound, z):
		
		size_k = 350 #size of the pt text file
		kpt = np.zeros((350))
		pt = np.zeros((350,self.lred))
		pt_temp = np.zeros((len(kbound),self.lred))
		#~ pt_final = np.zeros((len(kbound),len(np.atleast_1d(z))))
		pt_final = np.zeros((len(kbound)))
		for count,iz in enumerate(self.red):
			dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/'+pt_term+'_'+str(iz)+'.txt')

			with open(dat_file_path,'r') as f:  
				line = f.readline()
				for index_k in range(size_k):
					kpt[index_k] = float(line.split()[0])
					pt[index_k,count] = float(line.split()[1])
					line = f.readline()
			# interpolate on kbound
			pt_temp[:,count] = np.interp(kbound, kpt, pt[:,count]) 
			
			# interpolate on z
			for i in range(len(kbound)):
				#~ pt_final[i,:] = np.interp(z, self.red, pt_temp[i,:])
				pt_final[i] = np.interp(z, self.red, pt_temp[i,:])
				
		return pt_final
		
		
#-----------------------------------------------------------------------
		
	def bcoeff(self,z):
		b1 = np.zeros((self.lred))
		b2 = np.zeros((self.lred))
		b3 = np.zeros((self.lred))
		b4 = np.zeros((self.lred))
		b1_final = np.zeros(len(np.atleast_1d(z)))
		b2_final = np.zeros(len(np.atleast_1d(z)))
		b3_final = np.zeros(len(np.atleast_1d(z)))
		b4_final = np.zeros(len(np.atleast_1d(z)))
		
		if self.bmodel == 1:
			for count,iz in enumerate(self.red):
				dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/large_scale/'\
				'LS_z='+str(iz)+'_.txt')
				with open(dat_file_path,'r') as f: 
					line = f.readline()   
					b1[count] = float(line.split()[self.mbin])
					
			# interpolate on z
			b1_final = np.interp(z, self.red, b1)
			b2_final = np.interp(z, self.red, b2)
			b3_final = np.interp(z, self.red, b3)
			b4_final = np.interp(z, self.red, b4)
			
			return b1_final, b2_final, b3_final, b4_final # here b2_final, b3_final, b4_final == 0

		elif self.bmodel == 2:
			for count,iz in enumerate(self.red):
				dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV'\
				'/case'+str(self.kcase)+'/coeff_pl_0.0_z='+str(iz)+'.txt')
				with open(dat_file_path,'r') as f:
					for i, line in enumerate(f):
						if i == self.mbin: 
							b1[count] = float(line.split()[0])
							b2[count] = float(line.split()[1])
							b3[count] = float(line.split()[2])
							b4[count] = float(line.split()[3])
							
			# interpolate on z
			b1_final = np.interp(z, self.red, b1)
			b2_final = np.interp(z, self.red, b2)
			b3_final = np.interp(z, self.red, b3)
			b4_final = np.interp(z, self.red, b4)
					
			return b1_final, b2_final, b3_final, b4_final

		elif self.bmodel == 3:
			for count,iz in enumerate(self.red):
				dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV'\
				'/case'+str(self.kcase)+'/coeff_3exp_0.0_z='+str(iz)+'.txt')
				with open(dat_file_path,'r') as f: 
					for i, line in enumerate(f):
						if i == self.mbin: 
							b1[count] = float(line.split()[0])
							b2[count] = float(line.split()[1])
							b3[count] = float(line.split()[2])
							b4[count] = float(line.split()[3])
						
			# interpolate on z
			b1_final = np.interp(z, self.red, b1)
			b2_final = np.interp(z, self.red, b2)
			b3_final = np.interp(z, self.red, b3)
			b4_final = np.interp(z, self.red, b4)
					
			return b1_final, b2_final, b3_final, b4_final
		
#-----------------------------------------------------------------------

	def red_ps(self,data, kbound, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt, alpha):
			
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/TNS_coeff/'\
		'AB_'+str(self.mbin)+'_'+str(self.bmodel)+'_z='+str(self.z)+'.txt')	
		AB2 = np.zeros((len(kbound)))
		AB4 = np.zeros((len(kbound)))
		AB6 = np.zeros((len(kbound)))
		AB8 = np.zeros((len(kbound)))
		with open(dat_file_path,'r') as f: 
			line = f.readline()
			for index_k in range(len(kbound)):
				AB2[index_k] = float(line.split()[0])
				AB4[index_k] = float(line.split()[1])
				AB6[index_k] = float(line.split()[2])
				AB8[index_k] = float(line.split()[3])
				line = f.readline()
			
		if self.fog == 1:
			# compute the multipole expansion coefficients
			kappa = np.zeros((len(kbound)))
			coeffA = np.zeros((len(kbound)))
			coeffB = np.zeros((len(kbound)))
			coeffC = np.zeros((len(kbound)))
			coeffD = np.zeros((len(kbound)))
			coeffE = np.zeros((len(kbound)))
			
			sigma_v = (data.mcmc_parameters['sigma_v']['current'] *
			data.mcmc_parameters['sigma_v']['scale'])

			kappa = kbound*sigma_v*fz*Dz
			coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
			coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
			coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
			coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
			coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
	
			if self.rsd == 1:
				if self.bmodel == 1:
					b = b1 * alpha
					Pred = Pmod_dd*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
					
				elif self.bmodel == 2:
					b = (b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4))*alpha
					Pred = Pmod_dd*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
					
			elif self.rsd == 2:
				if self.bmodel == 1:
					b = b1 * alpha
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
					
				elif self.bmodel == 2:
					b = (b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4))*alpha
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
	
			elif self.rsd == 3:
				if self.bmodel == 1:
					b = b1
					Pred = b**2*Pmod_dd*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
					+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = b**2*Pmod_dd*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
					+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
					
				elif self.bmodel == 3:
					# here b3 == bs and b4 == b3nl
					PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
					2*b1*b4*F 
					PhhDT = b1* Pmod_dt + b2*G + b3*H + b4*F 
					Pred = PhhDD*coeffA  + 2/3.*fz*PhhDT*coeffB + 1/5.*fz**2*Pmod_tt*coeffC + 1/3.*AB2*coeffB \
					+ 1/5.*AB4*coeffC + 1/7.*AB6*coeffD + 1/9.*AB8*coeffE 
						
		#------------
		else:
			if self.rsd == 1:
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*(b**2 + 2/3.*b*fz + 1/5.*fz**2) 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = Pmod_dd*(b**2 + 2/3.*b*fz + 1/5.*fz**2 )
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
					
			elif self.rsd == 2:
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = Pmod_dd*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
	
			elif self.rsd == 3:
				if self.bmodel == 1:
					b = b1
					Pred = b**2*Pmod_dd + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
					+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = b**2*Pmod_dd + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
					+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
					
				elif self.bmodel == 3:
					# here b3 == bs and b4 == b3nl
					PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
					2*b1*b4*F 
					PhhDT = b1* Pmod_dt + b2*G + b3*H + b4*F 
					Pred = PhhDD  + 2/3.*fz*PhhDT + 1/5.*fz**2*Pmod_tt + 1/3.*AB2 \
					+ 1/5.*AB4 + 1/7.*AB6 + 1/9.*AB8 
			
		return Pred
#-----------------------------------------------------------------------

	def rescaling(self):
		
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/large_scale/'\
		'LS_z='+str(self.z)+'_.txt')
		with open(dat_file_path,'r') as f: 
				line = f.readline()   
				bcc_LS000 = float(line.split()[self.mbin])
				
		
		nu_masses = [0.03, 0.06, 0.1, 0.13, 0.15, 0.3, 0.45, 0.6]
		bcc_massive = np.zeros((len(nu_masses)))
		for count, mn in enumerate(nu_masses):
			dat_file_path2 = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/other neutrinos masses/'+\
			str(mn)+'eV/LS_z='+str(self.z)+'_.txt')
			with open(dat_file_path2,'r') as f2: 
				line2 = f2.readline()   
				bcc_massive[count] = float(line2.split()[self.mbin])
				
		# interpolate on z
		bcc_final = np.interp(self.Mnu, nu_masses, bcc_massive)/bcc_LS000

		return bcc_final
