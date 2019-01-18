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

		#~ #-------------------------------------------------
		#~ #---------------- Data ---------------------------
		#~ #-------------------------------------------------
		Class = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/class/test_z2_pk.dat')
		self.kclass = Class[:,0]
		self.Pclass = Class[:,1]
		
		
		# store the calibrated redshift
		self.red = [0.0,0.5,1.0,2.0] 
		self.lred = len(self.red) 
		
		
		#### load halo redshift space ps
		### for test
		self.data_directory = data.path['root']
		print self.data_directory
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/codes/analysis')
		sys.path.append(dat_file_path)
		from load_data import ld_data
		kcamb, Pcamb, k, Pmm, PH1, PH2, PH3 , PH4, errPhh1, errPhh2, errPhh3, errPhh4, bias1, bias2, bias3, bias4, \
		bias1s, bias2s, bias3s, bias4s, errb1, errb2, errb3, errb4, Pmono1, Pmono2, Pmono3, Pmono4, errPr1, errPr2, errPr3,\
		errPr4, kclass, Tm, Tcb, noise1, noise2, noise3, noise4 = ld_data(0.15, [0.0,0.5,1.0,2.0], 1)
		
		self.ksimu = k
		self.Psimu = Pmono1
		self.err = errPr1
		
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

		#### get the linear power spectrum from class
		#~ pk_lin = np.zeros((len(kclass)), 'float64')
		#~ for ik in xrange(len(kclass)):
			#~ pk_lin[ik] = cosmo.pk_lin(kclass[ik]*h, self.z)
		
		kbound = np.logspace(np.log10(self.kmin), np.log10(self.kmax), self.kbins)
		#### get the linear power spectrum from class. here multiply input k array by h because get_pk uses 1/mpc 
		pk_lin = cosmo.get_pk_array(kbound*h, redshift, len(kbound), znumber, 0) #if we want Pmm
		#~ pk_lin = cosmo.get_pk_cb_array(kclass, redshift, len(kclass), znumber, 0) # if we want Pcb
			
		
		### compare classy amplitude with classy
		import matplotlib.pyplot as plt
		#~ plt.plot(self.kclass, self.Pclass, c='r')
		#~ plt.plot(kbound, pk_lin*h**3, c='b')
		#~ plt.xscale('log')
		#~ plt.yscale('log')
		#~ plt.show()
		
		### rescale the amplitude of pk_lin accoridngly
		pk_lin *= h**3 
		
		#### get the non linear power spectrum from class
		#~ pk = np.zeros((len(kbound)), 'float64')
		#~ for ik in xrange(len(kbound)):
			#~ pk[ik] = cosmo.pk(kbound[ik], self.z)
				
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
		### compute the redshift power spectrum
		
		Pred = self.red_ps(data, kbound, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt)

		kill
			

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
		AB2, AB4, AB6, AB8 = fastpt.RSD_ABsum_components(pk_lin,fz,b1,C_window=C_window)

		
		
		#~ kill
		# compute ps in redshift space
		Pred = np.zeros((len(kbound)))
		Pred = PhhDD*coeffA  + 2/3.*fz*PhhDT*coeffB + 1/5.*fz**2*Pmod_tt*coeffC + 1/3.*AB2*coeffB \
		+ 1/5.*AB4*coeffC + 1/7.*AB6*coeffD + 1/9.*AB8*coeffE + A_shot
		
		
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
			
		print chi2
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
		pt_final = np.zeros((len(kbound),len(np.atleast_1d(z))))
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
				pt_final[i,:] = np.interp(z, self.red, pt_temp[i,:])
				
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
			print 'you chose the linear bias'
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
			print 'you chose the polynomial bias'
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
			print 'you chose the perturbation theory bias'
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

	def red_ps(self,data, kbound, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt):
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
				print 'you chose the Kaiser model'
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
					
			elif self.rsd == 2:
				print 'you chose the Scoccimaro model'
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = pk_lin*b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
	
			elif self.rsd == 3:
				print 'you chose the TNS model'
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						Pred[:,iz,count] = P_halo[:,iz,count]*coeffA[:,iz] + 2/3.*bmm[:,iz, count]*f[iz]*Pmod_dt[:,iz]*coeffB[:,iz] \
						+ 1/5.*f[iz]**2*Pmod_tt[:,iz]*coeffC[:,iz] + 1/3.*AB2[:,iz,count]*coeffB[:,iz] \
						+ 1/5.*AB4[:,iz,count]*coeffC[:,iz] + 1/7.*AB6[:,iz,count]*coeffD[:,iz] + 1/9.*AB8[:,iz,count]*coeffE[:,iz]
			elif self.rsd == 4:
				print 'you chose the eTNS model'
				dim = np.shape(P_halo)
				Pred = np.zeros(dim)
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						#~ ind2 = mbins.index(j)
						Pred[:,iz,count] = P_halo[:,iz,count]*coeffA[:,iz]  + 2/3.*f[iz]*PhhDT[:,iz,count]*coeffB[:,iz] \
						+ 1/5.*f[iz]**2*Pmod_tt[:,iz]*coeffC[:,iz] + 1/3.*AB2[:,iz,count]*coeffB[:,iz] \
						+ 1/5.*AB4[:,iz,count]*coeffC[:,iz] + 1/7.*AB6[:,iz,count]*coeffD[:,iz] + 1/9.*AB8[:,iz,count]*coeffE[:,iz]
		#------------
		else:
			if self.rsd == 1:
				print 'you chose the Kaiser model'
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						Pred[:,iz,count] = P_halo[:,iz,count] + 2/3.*bmm[:,iz, count]*f[iz] + 1/5.*f[iz]**2
			elif self.rsd == 2:
				print 'you chose the Scoccimaro model'
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						Pred[:,iz,count] = P_halo[:,iz,count] + 2/3.*bmm[:,iz, count]*f[iz]*Pmod_dt[:,iz]\
						+ 1/5.*f[iz]**2*Pmod_tt[:,iz]
			elif self.rsd == 3:
				print 'you chose the TNS model'
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						Pred[:,iz,count] = P_halo[:,iz,count]  + 2/3.*bmm[:,iz, count]*f[iz]*Pmod_dt[:,iz]\
						+ 1/5.*f[iz]**2*Pmod_tt[:,iz] + 1/3.*AB2[:,iz,count]+ 1/5.*AB4[:,iz,count]\
						+ 1/7.*AB6[:,iz,count]+ 1/9.*AB8[:,iz,count]
			elif self.rsd == 4:
				print 'you chose the eTNS model'
				dim = np.shape(P_halo)
				Pred = np.zeros(dim)
				for iz in xrange(znumber):
					for count,j in enumerate(Massbins):
						#~ ind2 = mbins.index(j)
						Pred[:,iz,count] = P_halo[:,iz,count]*coeffA[:,iz]  + 2/3.*f[iz]*PhhDT[:,iz,count]*coeffB[:,iz] \
						+ 1/5.*f[iz]**2*Pmod_tt[:,iz]*coeffC[:,iz] + 1/3.*AB2[:,iz,count]*coeffB[:,iz] \
						+ 1/5.*AB4[:,iz,count]*coeffC[:,iz] + 1/7.*AB6[:,iz,count]*coeffD[:,iz] + 1/9.*AB8[:,iz,count]*coeffE[:,iz]
			
		return Pred
	


