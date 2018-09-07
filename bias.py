
from classy import Class
from matplotlib.colors import LogNorm
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



def Halo(self, cosmo, data, model, case, Massbins):

	np.set_printoptions(precision=3)
	
	####################################################################
	#### check if the transfer functions were computed 
	test1 = data.cosmo_arguments
	output = test1.get('output')
	if 'mTk' not in output:
		raise ValueError('You forgot to declare mTk in Class output')

	####################################################################
	#### check if the non linear power spectrum has been requested 
        test2 = data.cosmo_arguments.keys()
	if 'non linear' not in test2:
		raise ValueError('You forgot to request the non linear spectrum')

	####################################################################
	#### check if kmax is defined in the data file of your likelihood 
	if 'z_max_pk' not in test2:
		raise ValueError('You must declare a z_max_pk for the computation of the transfer functions')

	####################################################################
	#### For now the simulations available are for Mv=0,0.06,0.10,0.15 
	#### Check if the total neutrino mass corresponds to one of the available ones
	#### get_current_derived_parameters returns a dict so must be converted
	m = cosmo.get_current_derived_parameters(['m_ncdm_tot'])
	m = m.values()
	m = [ round(elem, 2) for elem in m ]
	mv = [0.0, 0.03, 0.06, 0.10, 0.13, 0.15, 0.30]
	if m[0] not in mv:
		raise ValueError('Sorry the code is only available for Mv = 0.0, 0.03, 0.06, 0.10, 0.13, 0.15, 0.30 and your Mv is '+str(m[0])+'. Please modify you total neutrino mass.')

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
	#### Store the selected redshifts in a array and deduce its length for the loops
	#### array manipulation because len() and size only work for znumber >1
	a = np.array(redshift)
	znumber = a.size 
	redshift = np.zeros(znumber,'float64') 
	redshift[:] = a
	#~ print redshift

	#### CREATE FALSE REDSHIFT ARRAY FOR TEST
	redshift = [0.,0.5,1.,2.]
	znumber = len(redshift) # CHANGE WHEN THE CODE IS RUNNING
	print redshift
       
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
	
		b1 = np.zeros((len(red2),len(Massbins)))
		b2 = np.zeros((len(red2),len(Massbins)))
		bs = np.zeros((len(red2),len(Massbins)))
		b3nl = np.zeros((len(red2),len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_3exp_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					bs[ind,count] = f[ind2,2]
					b3nl[ind,count] = f[ind2,3]
					
	

	#-------------------------------------------------------------------
	if model =='pl':
	
		b1 = np.zeros((len(red2),len(Massbins)))
		b2 = np.zeros((len(red2),len(Massbins)))
		b3 = np.zeros((len(red2),len(Massbins)))
		b4 = np.zeros((len(red2),len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
				'/case2/coeff_pl_0.0_z='+str(i)+'.txt')
				f = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = f[ind2,0]
					b2[ind,count] = f[ind2,1]
					b3[ind,count] = f[ind2,2]
					b4[ind,count] = f[ind2,3]




	####################################################################
	#### and do a linear interpolation on requested redshifts
	#~ UNDEF = -999.0
	#~ b1_interp = np.interp(redshift, red1, b1, right=UNDEF)
	#~ b2_interp = np.interp(redshift, red1, b2, right=UNDEF)
	#~ b3_interp = np.interp(redshift, red1, b3, right=UNDEF)
	#~ b4_interp = np.interp(redshift, red1, b4, right=UNDEF)


	####################################################################
	#### get the rescaling coefficients according to neutrino mass
	bcc_LS000 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV/large_scale/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS000[ind,count] = f[ind2]
	#------------------------------
	bcc_LS003 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/other neutrinos masses/0.03/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS003[ind,count] = f[ind2]
	#------------------------------
	bcc_LS006 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/other neutrinos masses/0.06/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS006[ind,count] = f[ind2]
	#------------------------------
	bcc_LS010 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/other neutrinos masses/0.10/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS010[ind,count] = f[ind2]
	#------------------------------
	bcc_LS013 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/other neutrinos masses/0.13/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS013[ind,count] = f[ind2]
	#------------------------------
	bcc_LS015 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.15eV/large_scale/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS015[ind,count] = f[ind2]
	#------------------------------
	bcc_LS030 = np.zeros((len(red2),len(Massbins)))
	for i in red2:
		dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/other neutrinos masses/0.30/'\
		'LS_z='+str(i)+'_.txt')
		f = np.loadtxt(dat_file_path)
		ind = red2.index(i)
		for count,j in enumerate(Massbins):
			ind2 = mbins.index(j)
			bcc_LS030[ind,count] = f[ind2]

	
	#~ bcc_LS000 = np.interp(redshift, red2, bls000)
	#~ bcc_LS003 = np.interp(redshift, red2, bls003)
	#~ bcc_LS006 = np.interp(redshift, red2, bls006)
	#~ bcc_LS010 = np.interp(redshift, red2, bls010)
	#~ bcc_LS013 = np.interp(redshift, red2, bls013)
	#~ bcc_LS015 = np.interp(redshift, red2, bls015)
	#~ bcc_LS030 = np.interp(redshift, red2, bls030)

	
	####################################################################
	#### Since the cosmo.pk k's are bounded in [0.000000e+00:5.366287e+00]
	#### we must extract the k values from the get_transfer list so they coincide. 
	kget = cosmo.get_transfer(redshift[0])
	kclass = kget.get('k (h/Mpc)')
	bound = np.where((kclass > 0)&(kclass < 5.366287e+00))[0]
	kclass = kclass[bound]

	####################################################################
	#### get the index of the kmax value in the k.dot array and remove the elements out of range 
	#### if no kmax is defined the user is asked to define a kmax value in accordance with the 3 k-scale models in the article. 
	k_article = [0.35,0.2,0.15]

	max_z = np.max(redshift)
	UNDEF = -999.0
	if max_z > 3:
		kext = np.interp(max_z, red, [0.12,0.2,0.2,0.2], right=UNDEF)
	
	####################################################################
	#### select the k mode according to the ones in Raccanelli et al. 2017
	if case == 1:
		kmax = k_article[0]
	elif case == 2:
		kmax = k_article[1]
	elif case == 3:
		kmax = k_article[2]
		
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
	d_b = np.zeros((len(kclass), znumber), 'float64')
	d_cdm = np.zeros((len(kclass), znumber), 'float64')
	d_tot = np.zeros((len(kclass), znumber), 'float64')
	for i in xrange(znumber):
		transfer = cosmo.get_transfer(redshift[i])
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
	T_cb = np.zeros((len(kclass), znumber), 'float64')
	T_cb = (Omega_cdm * d_cdm + Omega_b * d_b)/(Omega_cdm + Omega_b)
    
    
	####################################################################
	#### get the non linear power spectrum from class
	pk = np.zeros((len(kclass), znumber), 'float64')
	for ik in xrange(len(kclass)):
		for iz in xrange(znumber):
			pk[ik,iz] = cosmo.pk(kclass[ik], redshift[iz])

	
	####################################################################
	#### get the linear power spectrum from class
	pk_lin = np.zeros((len(kclass), znumber), 'float64')
	for ik in xrange(len(kclass)):
		for iz in xrange(znumber):
			pk_lin[ik,iz] = cosmo.pk_lin(kclass[ik], redshift[iz])	
	
	####################################################################
	###### compute the one loop correction with FAST-PT for the expansion model


	if model == 'exp':
		#~ kclass = np.logspace(np.min(np.log10(kclasstemp)), np.max(np.log10(kclasstemp)), 122)
		Aprime = np.zeros((350,znumber))
		Bprime = np.zeros((350,znumber))
		Cprime = np.zeros((350,znumber))
		Dprime = np.zeros((350,znumber))
		Eprime = np.zeros((350,znumber))
		Fprime = np.zeros((350,znumber))
		Gprime = np.zeros((350,znumber))
		Hprime = np.zeros((350,znumber))
		Pmod_dd_prime = np.zeros((350,znumber))
		Pmod_dt_prime = np.zeros((350,znumber))
		Pmod_tt_prime = np.zeros((350,znumber))
		
		for count,iz in enumerate(redshift):
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/A_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Aprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/B_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Bprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/C_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Cprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/D_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Dprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/E_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Eprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/F_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Fprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/G_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Gprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/H_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Hprime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/Pmod_dd_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Pmod_dd_prime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/Pmod_dt_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Pmod_dt_prime[:,count] = f[:,1]
			#------------------------
			dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/0.0eV'\
			'/PT_coeff/Pmod_tt_'+str(iz)+'.txt')
			f = np.loadtxt(dat_file_path)
			kpt = f[:,0]
			Pmod_tt_prime[:,count] = f[:,1]
		
		
		### interpolate the pt coeff on the chosen scale
		A = np.zeros((len(kclass),znumber))
		B = np.zeros((len(kclass),znumber))
		C = np.zeros((len(kclass),znumber))
		D = np.zeros((len(kclass),znumber))
		E = np.zeros((len(kclass),znumber))
		F = np.zeros((len(kclass),znumber))
		G = np.zeros((len(kclass),znumber))
		H = np.zeros((len(kclass),znumber))
		Pmod_dd = np.zeros((len(kclass),znumber))
		Pmod_dt = np.zeros((len(kclass),znumber))
		Pmod_tt = np.zeros((len(kclass),znumber))
		for i in xrange(znumber):
			A[:,i] = np.interp(kclass, kpt, Aprime[:,i]) 
			B[:,i] = np.interp(kclass, kpt, Bprime[:,i]) 
			C[:,i] = np.interp(kclass, kpt, Cprime[:,i]) 
			D[:,i] = np.interp(kclass, kpt, Dprime[:,i]) 
			E[:,i] = np.interp(kclass, kpt, Eprime[:,i]) 
			F[:,i] = np.interp(kclass, kpt, Fprime[:,i]) 
			G[:,i] = np.interp(kclass, kpt, Gprime[:,i]) 
			H[:,i] = np.interp(kclass, kpt, Hprime[:,i]) 
			Pmod_dd[:,i] = np.interp(kclass, kpt, Pmod_dd_prime[:,i]) 
			Pmod_dt[:,i] = np.interp(kclass, kpt, Pmod_dt_prime[:,i]) 
			Pmod_tt[:,i] = np.interp(kclass, kpt, Pmod_tt_prime[:,i]) 
			
			
			
			
			
		#first mass range
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='+str(2.0)+'.txt')
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh1_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh2_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh3_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh4_realisation_0.15_z='+str(2.0)+'.txt')
		#~ k = d1[:,19]
		#~ Phh1 = np.zeros((len(k),10))
		#~ Phh2 = np.zeros((len(k),10))
		#~ Phh3 = np.zeros((len(k),10))
		#~ Phh4 = np.zeros((len(k),10))
		#~ Pshot1 = np.zeros((10))
		#~ Pshot2 = np.zeros((10))
		#~ Pshot3 = np.zeros((10))
		#~ Pshot4 = np.zeros((10))
		#~ pnum1 = [0,2,4,6,8,10,12,14,16,18]
		#~ pnum2 = [1,3,5,7,9,11,13,15,17,20]
		#~ for i in xrange(0,10):
			#~ Phh1[:,i]= d1[:,pnum1[i]]
			#~ Phh2[:,i]= d2[:,pnum1[i]]
			#~ Phh3[:,i]= d3[:,pnum1[i]]
			#~ Phh4[:,i]= d4[:,pnum1[i]]
			#~ Pshot1[i]= d1[0,pnum2[i]]
			#~ Pshot2[i]= d2[0,pnum2[i]]
			#~ Pshot3[i]= d3[0,pnum2[i]]
			#~ Pshot4[i]= d4[0,pnum2[i]]
			#~ Phh1[:,i] = Phh1[:,i]-Pshot1[i]
			#~ Phh2[:,i] = Phh2[:,i]-Pshot2[i]
			#~ Phh3[:,i] = Phh3[:,i]-Pshot3[i]
			#~ Phh4[:,i] = Phh4[:,i]-Pshot4[i]
			
		#~ PH1 = np.mean(Phh1[:,0:11], axis=1)
		#~ PH2 = np.mean(Phh2[:,0:11], axis=1)
		#~ PH3 = np.mean(Phh3[:,0:11], axis=1)
		#~ PH4 = np.mean(Phh4[:,0:11], axis=1)
		
		
		# compute the halo power spectrum given the coefficient
		PhhDD = np.zeros((len(kclass),znumber,len(Massbins)))
		PhhDT = np.zeros((len(kclass),znumber,len(Massbins)))
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				# density spectrum
				PhhDD[:,iz,count] = b1[iz,count]**2*Pmod_dd[:,iz] + b1[iz,count]*b2[iz,count]*A[:,iz] + 1/4.*b2[iz,count]**2*B[:,iz] + \
				b1[iz,count]*bs[iz,count]*C[:,iz] + 1/2.*b2[iz,count]*bs[iz,count]*D[:,iz] + 1/4.*bs[iz,count]**2*E[:,iz] +\
				2*b1[iz,count]*b3nl[iz,count]*F[:,iz] * (T_cb[:,iz]/d_tot[:,iz])**2
				# cross velocity spectrum
				PhhDT[:,iz,count] = b1[iz,count]* Pmod_dt[:,iz] + b2[iz,count]*G[:,iz] + bs[iz,count]*H[:,iz] + b3nl[iz,count] * F[:,iz]
					

		if m[0] == 0.03:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS003[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.06:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS006[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.10:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS010[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.13:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS013[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.15:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS015[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.30:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					PhhDD[:,iz,count] *= bcc_LS030[iz,count]/bcc_LS000[iz,count]

		
		
		
		# create a scale array limited by kmin and kmax
		try:
			self.kmax
		except:
			self.kmax = False
		if self.kmax:
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

		if self.kmin:
			lim_l = np.where(kclass >= self.kmin)[0]
		else:
			lim_l = np.where(kclass >= kmin)[0]
		#------------------------------------------------------------------
		##### =====> 
		kclass = kclass[lim_l[0]:klim_h+1]
		PhhDD = PhhDD[lim_l[0]:klim_h+1]
		PhhDT = PhhDT[lim_l[0]:klim_h+1]
		Pmod_tt = Pmod_tt[lim_l[0]:klim_h+1]

		#~ return kclass,PhhDD, PhhDT, Pmod_tt, k, PH1, PH2, PH3, PH4
		return kclass,PhhDD, PhhDT, Pmod_tt
		
		
	####################################################################
	###### compute the bias and halo power spectrum for power law model


	if model == 'pl':
		bcc = np.zeros((len(kclass), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bcc[:,iz, count] = b1[iz,count] + b2[iz,count]*(kclass**2) + b3[iz,count]*(kclass**3) \
				+ b4[iz,count]*(kclass**4) 
				
		if m[0] == 0.03:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS003[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.06:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS006[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.10:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS010[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.13:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS013[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.15:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS015[iz,count]/bcc_LS000[iz,count]
		elif m[0] == 0.30:
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					bcc[:,iz,count] *= bcc_LS030[iz,count]/bcc_LS000[iz,count]
		
				
		# compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.
		bmm = np.zeros((len(kclass), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bmm[:,iz, count] = bcc[:,iz,count] * (T_cb[:,iz]/d_tot[:,iz])# * (bcc_LS0[iz]/denom[iz])

		# Compute the halo Power spectrum in real space
		Phh = np.zeros((len(kclass), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
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
		if self.kmax:
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

		if self.kmin:
			lim_l = np.where(kclass >= self.kmin)[0]
		else:
			lim_l = np.where(kclass >= kmin)[0]
		#------------------------------------------------------------------
		##### =====> 
		kclass = kclass[lim_l[0]:klim_h+1]
		Phh = Phh[lim_l[0]:klim_h+1]
		
		#first mass range
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='+str(2.0)+'.txt')
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh1_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh2_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh3_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh4_realisation_0.15_z='+str(2.0)+'.txt')
		
		#~ k = d1[:,19]
		#~ Phh1 = np.zeros((len(k),10))
		#~ Phh2 = np.zeros((len(k),10))
		#~ Phh3 = np.zeros((len(k),10))
		#~ Phh4 = np.zeros((len(k),10))
		#~ Pshot1 = np.zeros((10))
		#~ Pshot2 = np.zeros((10))
		#~ Pshot3 = np.zeros((10))
		#~ Pshot4 = np.zeros((10))
		#~ pnum1 = [0,2,4,6,8,10,12,14,16,18]
		#~ pnum2 = [1,3,5,7,9,11,13,15,17,20]
		#~ for i in xrange(0,10):
			#~ Phh1[:,i]= d1[:,pnum1[i]]
			#~ Phh2[:,i]= d2[:,pnum1[i]]
			#~ Phh3[:,i]= d3[:,pnum1[i]]
			#~ Phh4[:,i]= d4[:,pnum1[i]]
			#~ Pshot1[i]= d1[0,pnum2[i]]
			#~ Pshot2[i]= d2[0,pnum2[i]]
			#~ Pshot3[i]= d3[0,pnum2[i]]
			#~ Pshot4[i]= d4[0,pnum2[i]]
			#~ Phh1[:,i] = Phh1[:,i]-Pshot1[i]
			#~ Phh2[:,i] = Phh2[:,i]-Pshot2[i]
			#~ Phh3[:,i] = Phh3[:,i]-Pshot3[i]
			#~ Phh4[:,i] = Phh4[:,i]-Pshot4[i]
			
		#~ PH1 = np.mean(Phh1[:,0:11], axis=1)
		#~ PH2 = np.mean(Phh2[:,0:11], axis=1)
		#~ PH3 = np.mean(Phh3[:,0:11], axis=1)
		#~ PH4 = np.mean(Phh4[:,0:11], axis=1)
		
		
		
		#~ return kclass, Phh, k, PH1, PH2, PH3, PH4
		
		
		return kclass, Phh
		





