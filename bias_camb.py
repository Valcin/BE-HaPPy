#### Link to Alvise article 1704
#### Link to VILLAESCUSA -NAVARRO article 1708
#### Link to Hernandez article 1608
#### Link to Linder article 1211

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
import numpy as np
import camb
import sys,os


################################## INPUT ######################################
# neutrino parameters
hierarchy = 'degenerate' #'degenerate', 'normal', 'inverted'
Mnu       = 0.00  #eV
Nnu       = 0  #number of massive neutrinos
Neff      = 3.046

# cosmological parameters
h       = 0.6711
Omega_c = 0.2685 - Mnu/(93.14*h**2)
Omega_b = 0.049
Omega_k = 0.0
tau     = None

# initial P(k) parameters
ns           = 0.9624
As           = 2.13e-9
pivot_scalar = 0.05
pivot_tensor = 0.05

# redshifts and k-range
redshifts    = [0.0, 0.5, 1, 2, 3, 99] 
kmax         = 10.0
k_per_logint = 10

# dz, relative difference dz/z to compute growths
dz = 0.01
###############################################################################

# create a new redshift list to compute growth rates
zs = []
for z in redshifts:
    dz_abs = (1.0+z)*dz
    if z==0.0:
        zs.append(z);  zs.append(z+dz_abs)
    else:
        zs.append(z-dz_abs);  zs.append(z);  zs.append(z+dz_abs)
z_list = redshifts;  redshifts = zs


Omega_cb = Omega_c + Omega_b

pars = camb.CAMBparams()

# set accuracy of the calculation
pars.set_accuracy(AccuracyBoost=5.0, lSampleBoost=5.0, 
                  lAccuracyBoost=5.0, HighAccuracyDefault=True, 
                  DoLateRadTruncation=True)

# set value of the cosmological parameters
pars.set_cosmology(H0=h*100.0, ombh2=Omega_b*h**2, omch2=Omega_c*h**2, 
                   mnu=Mnu, omk=Omega_k, 
                   neutrino_hierarchy=hierarchy, 
                   num_massive_neutrinos=Nnu,
                   nnu=Neff,
                   tau=tau)
                   
# set the value of the primordial power spectrum parameters
pars.InitPower.set_params(As=As, ns=ns, 
                          pivot_scalar=pivot_scalar, 
                          pivot_tensor=pivot_tensor)

# set redshifts, k-range and k-sampling
pars.set_matter_power(redshifts=redshifts, kmax=kmax, 
                      k_per_logint=k_per_logint)

# compute results
results = camb.get_results(pars)

# get raw matter power spectrum and transfer functions with strange k-binning
#k, zs, Pk = results.get_linear_matter_power_spectrum()
#Tk        = (results.get_matter_transfer_data()).transfer_data

# interpolate to get Pmm, Pcc...etc
k, zs, Pkmm = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                npoints=400, var1=7, var2=7, 
                                                have_power_spectra=True, 
                                                params=None)

k, zs, Pkcc = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                npoints=400, var1=2, var2=2, 
                                                have_power_spectra=True, 
                                                params=None)

k, zs, Pkbb = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                npoints=400, var1=3, var2=3, 
                                                have_power_spectra=True, 
                                                params=None)

k, zs, Pkcb = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                npoints=400, var1=2, var2=3, 
                                                have_power_spectra=True, 
                                                params=None)

Pkcb = (Omega_c**2*Pkcc + Omega_b**2*Pkbb +\
        2.0*Omega_b*Omega_c*Pkcb)/Omega_cb**2

k, zs, Pknn = results.get_matter_power_spectrum(minkh=2e-5, maxkh=kmax, 
                                                npoints=400, var1=6, var2=6, 
                                                have_power_spectra=True, 
                                                params=None)


print pars

# get sigma_8 and Hz in km/s/(kpc/h)
s8 = np.array(results.get_sigma8())
Hz = results.hubble_parameter(99.0)
print 'H(z=99)      = %.4f km/s/(kpc/h)'%(Hz/1e3/h)
print 'sigma_8(z=0) = %.4f'%s8[-1]


# do a loop over all redshifts
for i,z in enumerate(zs):

    fout1 = 'Pk_mm_z=%.3f.txt'%z
    fout2 = 'Pk_cc_z=%.3f.txt'%z
    fout3 = 'Pk_bb_z=%.3f.txt'%z
    fout4 = 'Pk_cb_z=%.3f.txt'%z
    fout5 = 'Pk_nn_z=%.3f.txt'%z

    np.savetxt(fout1,np.transpose([k,Pkmm[i,:]]))
    np.savetxt(fout2,np.transpose([k,Pkcc[i,:]]))
    np.savetxt(fout3,np.transpose([k,Pkbb[i,:]]))
    np.savetxt(fout4,np.transpose([k,Pkcb[i,:]]))
    np.savetxt(fout5,np.transpose([k,Pknn[i,:]]))


    #fout = 'Pk_trans_z=%.3f.txt'%z
    # notice that transfer functions have an inverted order:i=0 ==>z_max
    #np.savetxt(fout,np.transpose([Tk[0,:,i],Tk[1,:,i],Tk[2,:,i],Tk[3,:,i],
    #                               Tk[4,:,i],Tk[5,:,i],Tk[6,:,i]]))


# compute growth rates
for z in z_list:
    
    dz_abs = (1.0+z)*dz
    for suffix in ['mm','cb','nn']:

        fout = 'f%s_z=%.3f.txt'%(suffix,z)
        f2   = 'Pk_%s_z=%.3f.txt'%(suffix,z+dz_abs)

        if z==0.0:
            f1 = 'Pk_%s_z=%.3f.txt'%(suffix,z);  fac = 1.0
            
        else:
            f1 = 'Pk_%s_z=%.3f.txt'%(suffix,z-dz_abs);  fac = 2.0
            
        k1,Pk1 = np.loadtxt(f1,unpack=True)
        k2,Pk2 = np.loadtxt(f2,unpack=True)

        if np.any(k1!=k2):
            print 'Error!'; sys.exit()

        f = -0.5*(1.0+z)*np.log(Pk2/Pk1)/(fac*dz_abs)
        np.savetxt(fout,np.transpose([k1,f]))

        os.system('rm '+f2)
        if z!=0.0:  os.system('rm '+f1)


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
	mv = [0.0, 0.06, 0.10,0.15]
	if m[0] not in mv:
		raise ValueError('Sorry the code is only available for Mv = 0.0,0.06,0.10,0.15 and your Mv is '+str(m[0])+'. Please modify you total neutrino mass.')

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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case1/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_3exp_'+str(m[0])+'_z='+str(i)+'.txt')
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case1/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
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
				dat_file_path = os.path.join(self.data_directory, 'montepython/BE_HaPPy/coefficients/'+str(m[0])+\
				'eV/case2/coeff_pl_'+str(m[0])+'_z='+str(i)+'.txt')
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

	#~ bcc_LS0 = np.interp(redshift, red2, bls0)
	#~ bcc_LS006 = np.interp(redshift, red2, bls006)
	#~ bcc_LS010 = np.interp(redshift, red2, bls010)
	#~ bcc_LS015 = np.interp(redshift, red2, bls015)

	
	####################################################################
	#### Since the cosmo.pk k's are bounded in [0.000000e+00:5.366287e+00]
	#### we must extract the k values from the get_transfer list so they coincide. 
	kget = cosmo.get_transfer(redshift[0])
	kclasstemp = kget.get('k (h/Mpc)')
	bound = np.where((kclasstemp > 0)&(kclasstemp < 5.366287e+00))[0]
	kclasstemp = kclasstemp[bound]

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
	d_b = np.zeros((len(kclasstemp), znumber), 'float64')
	d_cdm = np.zeros((len(kclasstemp), znumber), 'float64')
	d_tot = np.zeros((len(kclasstemp), znumber), 'float64')
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
	T_cb = np.zeros((len(kclasstemp), znumber), 'float64')
	T_cb = (Omega_cdm * d_cdm + Omega_b * d_b)/(Omega_cdm + Omega_b)
    
    
	####################################################################
	#### get the non linear power spectrum from class and rescale it in (Mpc/h)**3
	pk = np.zeros((len(kclasstemp), znumber), 'float64')
	for ik in xrange(len(kclasstemp)):
		for iz in xrange(znumber):
			pk[ik,iz] = cosmo.pk(kclasstemp[ik], redshift[iz])

	
	####################################################################
	#### get the linear power spectrum from class
	pk_lin = np.zeros((len(kclasstemp), znumber), 'float64')
	for ik in xrange(len(kclasstemp)):
		for iz in xrange(znumber):
			pk_lin[ik,iz] = cosmo.pk_lin(kclasstemp[ik], redshift[iz])	
	
	####################################################################
	###### compute the one loop correction with FAST-PT for the expansion model


	if model == 'exp':
		# evenly logged kclass and interpolate pk_lin for fast-PT
		kclass = np.logspace(np.min(np.log10(kclasstemp)), np.max(np.log10(kclasstemp)), 122)
		pk_lin2 = np.zeros((len(kclass),znumber))
		T_cb2 = np.zeros((len(kclass),znumber))
		d_tot2 = np.zeros((len(kclass),znumber))
		
		for iz in xrange(znumber):

			# set the parameters for the power spectrum window and
			# Fourier coefficient window 
			C_window=.75

			# interpolate on the new scale array
			pk_lin2[:,iz] = np.interp(kclass,kclasstemp, pk_lin[:,iz])
			T_cb2[:,iz] = np.interp(kclass,kclasstemp, T_cb[:,iz])
			d_tot2[:,iz] = np.interp(kclass,kclasstemp, d_tot[:,iz])
			
			# padding length 
			nu=-2; n_pad=len(kclass)
			n_pad=int(0.5*len(kclass))
			to_do=['all']
							
			
			
			# initialize the FASTPT class 
			# including extrapolation to higher and lower k  
			# time the operation
			fastpt=FPT.FASTPT(kclass,to_do=to_do,n_pad=n_pad) 
				
			# calculate 1loop SPT (and time the operation) for density
			P_spt_dd=fastpt.one_loop_dd(pk_lin2[:,iz],C_window=C_window)
				
			
			# calculate 1loop SPT (and time the operation) for velocity
			P_spt_tt=fastpt.one_loop_tt(pk_lin2[:,iz],C_window=C_window)

				
			# calculate 1loop SPT (and time the operation) for velocity - density
			P_spt_dt=fastpt.one_loop_dt(pk_lin2[:,iz],C_window=C_window)

				
			#calculate tidal torque EE and BB P(k)
			#~ P_RSD=fastpt.RSD_components(P,1.0,C_window=C_window)	

			Pmod_tt = np.zeros((len(kclass),znumber))
			# update the power spectrum
			Pmod_dd=pk_lin2[:,iz]+P_spt_dd[0]
			Pmod_dt=pk_lin2[:,iz]+P_spt_dt[0]
			Pmod_tt[:,iz]=pk_lin2[:,iz]+P_spt_tt[0]
			
			A = P_spt_dd[2]
			B = P_spt_dd[3]
			C = P_spt_dd[4]
			D = P_spt_dd[5]
			E = P_spt_dd[6]
			F = P_spt_dd[7]
			G = P_spt_dt[2]
			H = P_spt_dt[3]
			
			
			
			#first mass range
			#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='+str(2.0)+'.txt')
			#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='+str(2.0)+'.txt')
			#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='+str(2.0)+'.txt')
			#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='+str(2.0)+'.txt')
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
			
			
			# compute the halo power spectrum given the coefficient
			PhhDD = np.zeros((len(kclass),znumber,len(Massbins)))
			PhhDT = np.zeros((len(kclass),znumber,len(Massbins)))
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					# density spectrum
					PhhDD[:,iz,count] = b1[iz,count]**2*Pmod_dd + b1[iz,count]*b2[iz,count]*A + 1/4.*b2[iz,count]**2*B + \
					b1[iz,count]*bs[iz,count]*C + 1/2.*b2[iz,count]*bs[iz,count]*D + 1/4.*bs[iz,count]**2*E +\
					2*b1[iz,count]*b3nl[iz,count]*F * (T_cb2[:,iz]/d_tot2[:,iz])**2
					# cross velocity spectrum
					PhhDT[:,iz,count] = b1[iz,count]* Pmod_dt + b2[iz,count]*G + bs[iz,count]*H + b3nl[iz,count] * F
					
		
		
		# rescale the k and power spectrum because of classy/class difference (explain this power of k ?)
		kclass /= h
		# For z == 0.0:
		PhhDD[:,0,:] /= h**(3/2.)
		PhhDT[:,0,:] /= h**(3/2.)
		Pmod_tt[:,0] /= h**(3/2.)
		# For z == 0.5:
		PhhDD[:,1,:] /= 1.
		PhhDT[:,1,:] /= 1.
		Pmod_tt[:,1] /= 1.
		# For z == 1.0:
		PhhDD[:,2,:] *= h
		PhhDT[:,2,:] *= h
		Pmod_tt[:,2] *= h
		# For z == 2.0:
		PhhDD[:,3,:] *= h**3
		PhhDT[:,3] *= h**3
		#~ Pmod_tt[:,3,:] *= h**3

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

		#~ print np.max(kclasstemp)
		return kclass,PhhDD, PhhDT, Pmod_tt
		
		
	####################################################################
	###### compute the bias and halo power spectrum for power law model


	if model == 'pl':
		bcc = np.zeros((len(kclasstemp), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bcc[:,iz, count] = b1[iz,count] + b2[iz,count]*(kclasstemp**2) + b3[iz,count]*(kclasstemp**3) \
				+ b4[iz,count]*(kclasstemp**4) 
				
		# compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.
		bmm = np.zeros((len(kclasstemp), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				bmm[:,iz, count] = bcc[:,iz,count] * (T_cb[:,iz]/d_tot[:,iz])# * (bcc_LS0[iz]/denom[iz])

		# Compute the halo Power spectrum in real space
		Phh = np.zeros((len(kclasstemp), znumber, len(Massbins)), 'float64')
		for iz in xrange(znumber):
			for count,j in enumerate(Massbins):
				ind2 = mbins.index(j)
				Phh[:,iz,count] = pk[:,iz] * bmm[:,iz, count]**2
			
			
		# rescale the k and power spectrum because of classy/class difference
		kclasstemp /= h
		Phh *= h**3
		
		# create a scale array limited by kmin and kmax
		try:
			self.kmax
		except:
			self.kmax = False
		if self.kmax:
			lim_h = np.where(kclasstemp <= self.kmax)[0]
		else:
			lim_h = np.where(kclasstemp <= kmax)[0]

		klim_h = np.amax(lim_h) # define the higher limit
		
		#-------------------------------------------------------------------
		Vs = 1000**3 # for a periodic box of 1000 h-1Mpc
		kmin = 2 * math.pi * Vs**(-1/3.)

		try:
			self.kmin
		except:
			self.kmin = False

		if self.kmin:
			lim_l = np.where(kclasstemp >= self.kmin)[0]
		else:
			lim_l = np.where(kclasstemp >= kmin)[0]
		#------------------------------------------------------------------
		##### =====> 
		kclasstemp = kclasstemp[lim_l[0]:klim_h+1]
		Phh = Phh[lim_l[0]:klim_h+1]
		
		return kclasstemp, Phh
		
		

	##### check the total neutrino mass and get the denom accordingly
	#~ if m[0] == 0.:
		#~ print 'Mv = '+str(m[0])
		#~ denom = bcc_LS0
	#~ elif m[0] == 0.06:
		#~ print 'Mv = '+str(m[0])
		#~ denom = bcc_LS006
	#~ elif m[0] == 0.10:
		#~ print 'Mv = '+str(m[0])
		#~ denom = bcc_LS010
	#~ elif m[0] == 0.15:
		#~ print 'Mv = '+str(m[0])
		#~ denom = bcc_LS015


	

	

	



