from classy import Class
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import scipy.constants as const
import math
import numpy as np




def Finger(cosmo,redshift,znumber, k_dot):

	Vs = 1000**3 # for a periodic box of 1000 h-1Mpc
	
	#### import Omega_b and Omega_cdm from class. Remember to add Omega_cdm in classy and recompile after
	Omega_b = cosmo.Omega_b()
	Omega_cdm = cosmo.Omega_cdm()
	Omega_m = cosmo.Omega_m()
	
	#### get the value of h to rescale the power spectrum and wavenumber
	param = cosmo.get_current_derived_parameters(['h','Omega0_lambda'])
	param = param.values()
	param = [ round(elem, 5) for elem in param ]
	h = param[0]
	Omega_lambda = param[1]

	#### Define the linear growth factor and growth rate (growth factor f in class)
	f = np.zeros(znumber, 'float64')
	D = np.zeros(znumber, 'float64')
	for i in xrange(znumber):
		f[i] = cosmo.scale_independent_growth_factor_f(redshift[i])
		D[i] = cosmo.scale_independent_growth_factor(redshift[i])

	#### Query the cosmo module for the Hubble rate (in 1/Mpc) and convert it to km/s/Mpc
	#### compute a velocity dispersion to simulate F.o.G effects
	h_z = np.zeros(znumber, 'float64')
	v_disp = np.zeros(znumber, 'float64')
	FoG = np.zeros((len(k_dot), znumber), 'float64')
	c_light_km_per_sec = const.c/1000.
	G = 4.302e-9 # in Mpc Msol-1 (km/s)**2
	from hmf import MassFunction
	hmf = MassFunction(Mmin = 13.5, lnk_min = -2.22, lnk_max = -0.30)

	
	for i in xrange(znumber):
		hmf.update(z=redshift[i])
		h_z[i] = cosmo.Hubble(redshift[i]) * c_light_km_per_sec / 100
		Omeg_m_z = Omega_m * (1 + redshift[i])**3 / (Omega_m * (1 + redshift[i])**3 + Omega_lambda) # with Omega_k = 0
		mass_func = hmf.dndm
		mass = hmf.m
		rad = hmf.radii
		#~ plt.plot(ma,mass_func, label=str(redshift[i]))
		n = mass * mass_func* Vs
		b = np.average(mass, weights=n)
		r_vir = (0.784 *(b/ 1e8)**(1/3.) * (Omega_m / Omeg_m_z * 200 / 18 / math.pi**2)**(-1/3.)* ((1+redshift[i])/10)**(-1))/1000 #in mpc
		M_200 = 100 * r_vir**3 * (h_z[i]*100)**2 / G	#/ (1+redshift[i])**3
		v_disp[i] = 10**(np.log10(1082.9) + (0.3361)* np.log10(h_z[i]* M_200 / 1e15))
		#~ v_disp[i] = 1019*(h_z[i]* M_200 / 1e15)**(1/3.)
		
	
		v_disp[i] = v_disp[i]/(h*100) # from km/s to Mpc/h
		FoG[:,i] = np.exp(-(k_dot* v_disp[i]*f[i]/h_z[i]**3)**2)
		
	
	#~ print v_disp
	
	return FoG
	
####------------------------------------------------------------------------------------
	
def Lin_Kaiser(cosmo,redshift,znumber, ksize, bias):
		
	#### Define the linear growth rate (growth factor f in class)
	f = np.zeros(znumber, 'float64')
	mono = np.zeros((ksize, znumber), 'float64')
	quadru = np.zeros((ksize, znumber), 'float64')
	beta = np.zeros((ksize, znumber), 'float64')
		 
	for i in xrange(znumber):
		f[i] = cosmo.scale_independent_growth_factor_f(redshift[i])
		if bias is None:
			mono[:,i] = (1 + 2/3.*(f[i]) + 1/5.*(f[i])**2) 
			quadru[:,i] = (4/3.*(f[i]) + 4/7.*(f[i])**2)
		else:
			beta[:,i] = f[i]/bias[:,i]
			mono[:,i] = (1 + 2/3.*(beta[:,i]) + 1/5.*(beta[:,i])**2) 
			quadru[:,i] = (4/3.*(beta[:,i]) + 4/7.*(beta[:,i])**2)	
	return mono, quadru
	
####--------------------------------------------------------------------------------------

	
def NonLin_Kaiser(cosmo,redshift,znumber, ksize, bias):
	
	
	
	
	
	
	return


####---------------------------------------------------------------------------------------
	
def TNS(cosmo,redshift,znumber, ksize, bias):
	
	
	
	
	return TNS
