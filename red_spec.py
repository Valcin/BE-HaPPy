#### Link to Alvise article 1704
#### Link to VILLAESCUSA -NAVARRO article 1708
#### Link to Hernandez article 1608
#### Link to Linder article 1211

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
from bias import Halo

def rspec(self, cosmo, data, model, case, Massbins, RSD=None, fog=None):

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
	
	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)

	####################################################################
	#### Store the mass bins available in an array
	mbins = ['M1', 'M2', 'M3', 'M4']

	####################################################################
	#### Define the linear growth factor and growth rate (growth factor f in class)
	f = np.zeros(znumber, 'float64')
	D = np.zeros(znumber, 'float64')
	for iz in xrange(znumber):
		f[iz] = cosmo.scale_independent_growth_factor_f(redshift[iz])
		D[iz] = cosmo.scale_independent_growth_factor(redshift[iz])
		
	print f, D
		
	####################################################################
	#### get the value of h to rescale the power spectrum and wavenumber
	#### because of the difference between Class and Class python wrapper 
	param = cosmo.get_current_derived_parameters(['h','Omega0_lambda'])
	param = param.values()
	param = [ round(elem, 5) for elem in param ]
	h = param[0]
	Omega_lambda = param[1]

	####################################################################
	#### For now the simulations available are for Mv=0,0.06,0.10,0.15 
	#### Check if the total neutrino mass corresponds to one of the available ones
	#### get_current_derived_parameters returns a dict so must be converted
	m = cosmo.get_current_derived_parameters(['m_ncdm_tot'])
	m = m.values()
	m = [ round(elem, 2) for elem in m ]
	mv = [0.0, 0.06, 0.10,0.15]
	if m[0] not in mv:
		raise ValueError('Sorry the code is only available for Mv = 0.0,0.06,0.10,0.15 and your Mv is '+str(m[0])+'. Please modify you total neutrino mass.')

	
	####################################################################
	#### get the coefficients from the dat files in the data directory 
	#### get the large scale amplitude of bcc at different z and for different neutrino masses
	self.data_directory = data.path['root']

	#-------------------------------------------------------------------
	if model =='exp':
	
		b1 = np.zeros((len(red2),len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case1/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]


		if case == 2:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]


		if case == 3:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]

	#-------------------------------------------------------------------
	if model =='pl':
	
		b1 = np.zeros((len(red2),len(Massbins)))

		if case == 1:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case1/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]


		if case == 2:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]


		if case == 3:
			for i in red2:
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
				fichier = np.loadtxt(dat_file_path)
				ind = red2.index(i)
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					b1[ind,count] = fichier[ind2,0]


	####################################################################
	#### compute the fog multipole coefficient if requested
	if fog:
		kappa = k*(data.mcmc_parameters['sigma_v']['current']*data.mcmc_parameters['sigma_v']['scale'])
		coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
	
	####################################################################
	if model == 'exp':
		k, Pdd,Pdt,Ptt = Halo(self, cosmo, data, model, case, Massbins)
		
		# get the linear power spectrum from class
		pk_lin = np.zeros((len(k), znumber), 'float64')
		for ik in xrange(len(k)):
			for iz in xrange(znumber):
				pk_lin[ik,iz] = cosmo.pk_lin(k[ik], redshift[iz])	
		
		# set the parameters for the power spectrum window and
		# Fourier coefficient window 
		C_window=.75

		# padding length 
		nu=-2; n_pad=int(0.5*len(k))
		to_do=['all']
		
		# evenly logged kclass and interpolate pk_lin for fast-PT
		kbis = np.logspace(np.min(np.log10(k)), np.max(np.log10(k)), 122)
		pk_lin2 = np.zeros((len(kbis),znumber))
		Pdd2 = np.zeros((len(kbis),znumber,len(Massbins)))
		Pdt2 = np.zeros((len(kbis),znumber,len(Massbins)))
		Ptt2 = np.zeros((len(kbis),znumber))
		


		# initialize the FASTPT class 
		# including extrapolation to higher and lower k  
		# time the operation
		fastpt=FPT.FASTPT(kbis,to_do=to_do,n_pad=n_pad) 
		# taruya multipole coefficient
		
		
		
		if fog:
			# compute the halo ps in redshift space
			Phh = np.zeros((len(kbis),znumber,len(Massbins)))
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					# interpolate on the new scale array
					pk_lin2[:,iz] = np.interp(kbis,k, pk_lin[:,iz])
					Pdd2[:,iz,count] = np.interp(kbis,k, Pdd[:,iz,count])
					Pdt2[:,iz,count] = np.interp(kbis,k, Pdt[:,iz,count])
					Ptt2[:,iz] = np.interp(kbis,k, Ptt[:,iz])
					AB2,AB4,AB6,AB8 = fastpt2.RSD_ABsum_components(pk_lin2[:,iz],f[iz], b1[iz,count],C_window=C_window) #tns coeff
					Phh[:,iz,count] = Pdd2[:,iz,count]*coeffA + 2/3.*f[iz]*Pdt2[:,iz,count]*coeffB +\
					1/5.*f[iz]**2*Ptt2[:,iz]*coeffC + 1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE
		else:
			Phh = np.zeros((len(kbis),znumber,len(Massbins)))
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					# interpolate on the new scale array
					pk_lin2[:,iz] = np.interp(kbis,k, pk_lin[:,iz])
					Pdd2[:,iz,count] = np.interp(kbis,k, Pdd[:,iz,count])
					Pdt2[:,iz,count] = np.interp(kbis,k, Pdt[:,iz,count])
					Ptt2[:,iz] = np.interp(kbis,k, Ptt[:,iz])
					AB2,AB4,AB6,AB8 = fastpt.RSD_ABsum_components(pk_lin2[:,iz],f[iz], b1[iz,count],C_window=C_window) #tns coeff
					Phh[:,iz,count] = Pdd2[:,iz,count] + 2/3.*f[iz]*Pdt2[:,iz,count] +\
					1/5.*f[iz]**2*Ptt2[:,iz] + 1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6 + 1/9.*AB8
					
		
		#---------------------------------------------------
		#--------- halo redshift space -------------------------
		#---------------------------------------------------
		d1a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_0_z='+str(2.0)+'.txt')
		d1b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_0_z='+str(2.0)+'.txt')
		d1c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_0_z='+str(2.0)+'.txt')
		d1d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_0_z='+str(2.0)+'.txt')
		d2a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_1_z='+str(2.0)+'.txt')
		d2b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_1_z='+str(2.0)+'.txt')
		d2c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_1_z='+str(2.0)+'.txt')
		d2d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_1_z='+str(2.0)+'.txt')
		d3a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_2_z='+str(2.0)+'.txt')
		d3b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_2_z='+str(2.0)+'.txt')
		d3c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_2_z='+str(2.0)+'.txt')
		d3d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_2_z='+str(2.0)+'.txt')


		kx1a = d1a[:,19]
		Px1a = np.zeros((len(kx1a),10))
		Px1b = np.zeros((len(kx1a),10))
		Px1c = np.zeros((len(kx1a),10))
		Px1d = np.zeros((len(kx1a),10))
		Pxshot1a = np.zeros((10))
		Pxshot1b = np.zeros((10))
		Pxshot1c = np.zeros((10))
		Pxshot1d = np.zeros((10))
		pnum1 = [0,2,4,6,8,10,12,14,16,18]
		pnum2 = [1,3,5,7,9,11,13,15,17,20]
		for i in xrange(0,10):
			Px1a[:,i]= d1a[:,pnum1[i]]
			Px1b[:,i]= d1b[:,pnum1[i]]
			Px1c[:,i]= d1c[:,pnum1[i]]
			Px1d[:,i]= d1d[:,pnum1[i]]
			Pxshot1a[i]= d1a[0,pnum2[i]]
			Pxshot1b[i]= d1b[0,pnum2[i]]
			Pxshot1c[i]= d1c[0,pnum2[i]]
			Pxshot1d[i]= d1d[0,pnum2[i]]
		kx2a = d2a[:,19]
		Px2a = np.zeros((len(kx2a),10))
		Px2b = np.zeros((len(kx2a),10))
		Px2c = np.zeros((len(kx2a),10))
		Px2d = np.zeros((len(kx2a),10))
		Pxshot2a = np.zeros((10))
		Pxshot2b = np.zeros((10))
		Pxshot2c = np.zeros((10))
		Pxshot2d = np.zeros((10))
		pnum1 = [0,2,4,6,8,10,12,14,16,18]
		pnum2 = [1,3,5,7,9,11,13,15,17,20]
		for i in xrange(0,10):
			Px2a[:,i]= d2a[:,pnum1[i]]
			Px2b[:,i]= d2b[:,pnum1[i]]
			Px2c[:,i]= d2c[:,pnum1[i]]
			Px2d[:,i]= d2d[:,pnum1[i]]
			Pxshot2a[i]= d2a[0,pnum2[i]]
			Pxshot2b[i]= d2b[0,pnum2[i]]
			Pxshot2c[i]= d2c[0,pnum2[i]]
			Pxshot2d[i]= d2d[0,pnum2[i]]
		kx3a = d3a[:,19]
		Px3a = np.zeros((len(kx3a),10))
		Px3b = np.zeros((len(kx3a),10))
		Px3c = np.zeros((len(kx3a),10))
		Px3d = np.zeros((len(kx3a),10))
		Pxshot3a = np.zeros((10))
		Pxshot3b = np.zeros((10))
		Pxshot3c = np.zeros((10))
		Pxshot3d = np.zeros((10))
		pnum1 = [0,2,4,6,8,10,12,14,16,18]
		pnum2 = [1,3,5,7,9,11,13,15,17,20]
		for i in xrange(0,10):
			Px3a[:,i]= d3a[:,pnum1[i]]
			Px3b[:,i]= d3b[:,pnum1[i]]
			Px3c[:,i]= d3c[:,pnum1[i]]
			Px3d[:,i]= d3d[:,pnum1[i]]
			Pxshot3a[i]= d3a[0,pnum2[i]]
			Pxshot3b[i]= d3b[0,pnum2[i]]
			Pxshot3c[i]= d3c[0,pnum2[i]]
			Pxshot3d[i]= d3d[0,pnum2[i]]
			
		for i in xrange(0,10):
			Px1a[:,i] = Px1a[:,i]-Pxshot1a[i]
			Px1b[:,i] = Px1b[:,i]-Pxshot1b[i]
			Px1c[:,i] = Px1c[:,i]-Pxshot1c[i]
			Px1d[:,i] = Px1d[:,i]-Pxshot1d[i]
			Px2a[:,i] = Px2a[:,i]-Pxshot2a[i]
			Px2b[:,i] = Px2b[:,i]-Pxshot2b[i]
			Px2c[:,i] = Px2c[:,i]-Pxshot2c[i]
			Px2d[:,i] = Px2d[:,i]-Pxshot2d[i]
			Px3a[:,i] = Px3a[:,i]-Pxshot3a[i]
			Px3b[:,i] = Px3b[:,i]-Pxshot3b[i]
			Px3c[:,i] = Px3c[:,i]-Pxshot3c[i]
			Px3d[:,i] = Px3d[:,i]-Pxshot3d[i]
			
			nul1a = np.where(Px1a[:,i] < 0)[0]
			Px1a[nul1a,i] = 0
			nul1b = np.where(Px1b[:,i] < 0)[0]
			Px1b[nul1b,i] = 0
			nul1c = np.where(Px1c[:,i] < 0)[0]
			Px1c[nul1c,i] = 0
			nul1d = np.where(Px1d[:,i] < 0)[0]
			Px1d[nul1d,i] = 0
			nul2a = np.where(Px2a[:,i] < 0)[0]
			Px2a[nul2a,i] = 0
			nul2b = np.where(Px2b[:,i] < 0)[0]
			Px2b[nul2b,i] = 0
			nul2c = np.where(Px2c[:,i] < 0)[0]
			Px2c[nul2c,i] = 0
			nul2d = np.where(Px2d[:,i] < 0)[0]
			Px2d[nul2d,i] = 0
			nul3a = np.where(Px3a[:,i] < 0)[0]
			Px3a[nul3a,i] = 0
			nul3b = np.where(Px3b[:,i] < 0)[0]
			Px3b[nul3b,i] = 0
			nul3c = np.where(Px3c[:,i] < 0)[0]
			Px3c[nul3c,i] = 0
			nul3d = np.where(Px3d[:,i] < 0)[0]
			Px3d[nul3d,i] = 0

		Pmono1temp = (Px1a + Px2a + Px3a)/3
		Pmono2temp = (Px1b + Px2b + Px3b)/3
		Pmono3temp = (Px1c + Px2c + Px3c)/3
		Pmono4temp = (Px1d + Px2d + Px3d)/3


		### do the mean and std over quantitites ###
		
		Pmono1 = np.mean(Pmono1temp[:,0:11], axis=1)
		Pmono2 = np.mean(Pmono2temp[:,0:11], axis=1)
		Pmono3 = np.mean(Pmono3temp[:,0:11], axis=1)
		Pmono4 = np.mean(Pmono4temp[:,0:11], axis=1)

		return kbis,Phh, kx1a, Pmono1, Pmono2, Pmono3, Pmono4

	



