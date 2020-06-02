### code written by David Valcin
### for calibration details see arXiv:1901.06045

from classy import Class
import numpy as np
import sys
sys.path.append('/home/david/codes/BE_HaPPy') # the folder where BE-Happy is installed
from main import ps_calc
import matplotlib.pyplot as plt

########################################################################
########################################################################
### BE_HaPPY configuration
########################################################################
########################################################################

### if you want to compare several spectra define quantities as list (e.g. coord = [0, 1]) and compute for the desired index

# real space or redshift space(0 or 1), 1 for redshift space
coord = [0,1]
#choice of k model between (1,2,3) cf paper
kcase = '2'
# integration boundaries for k (in h/Mpc).
# It has to be within the boundaries of the k model above, respectively 0.15, 0.2, 0.4
kmin = 0.0001
kmax = 0.2
kbins = 261

# desired redshift must be between [0-2]
z = 1.5
# choice of neutrino mass between [0.0 - 0.6]eV
Mnu = 0.0
#choice of mass bins between (0,1,2,3) cf paper
mbin = 0
# choice of bias model between (1,2,3) cf paper: linear, polynomial, perturbation theory
bias_model = 3
#choice of RSD model between (1,2,3) cf paper
rsd = 3
# whether or not to include FoG (0 or 1), 1 for yes. sigma_v is defined with a default value fitted with the
# coefficients but the user can use different values
fog = 0
#sigma_v = [7.0, 7.0, 7.0, 7.0]



########################################################################
########################################################################
### Here an example
### Import a linear power spectrum. Make sure to use values close the paper reference cosmology
########################################################################
########################################################################
znumber = 1
zmax = 2 # needed for transfer

params = {'output': 'mPk mTk', 'z_max_pk': zmax, 'non linear': 'halofit', 'sigma8': 0.79,
		'n_s': 0.95, 'h': 0.7, 'omega_b': 0.023, 'Omega_cdm': 0.1093, 'Omega_Lambda' : 0.73}

# Create an instance of the CLASS wrapper
cosmo = Class()

# Set the parameters to the cosmological code
cosmo.set(params)
cosmo.compute()

#### Define the linear growth factor and growth rate (growth factor f in class)
h = cosmo.h()
fz = cosmo.scale_independent_growth_factor_f(z)
Dz = cosmo.scale_independent_growth_factor(z)
karray = np.logspace(np.log10(kmin), np.log10(kmax), kbins)
#~ Plin = cosmo.get_pk_cb_array(karray*h, np.array([z]), len(karray), znumber, 0) # if we want Pcb
Plin = np.zeros(len(karray))
for i,k in enumerate(karray) :
	Plin[i] = (cosmo.pk(k ,z)) # function .pk(k,z)
	
########################################################################
########################################################################
### Compute the non linear power spectrum
########################################################################
########################################################################



# A summary of your configuration is printed in the terminal
P1 = ps_calc(coord[0],kcase, Mnu, mbin, rsd, bias_model, karray, z, fog, Plin, fz, Dz)
P2 = ps_calc(coord[1],kcase, Mnu, mbin, rsd, bias_model, karray, z, fog, Plin, fz, Dz)


### save results
#~ for j in range(len(karray)):
	#~ with open('Pk_z'+str(z)+'_M'+str(mbin+1)+'.txt', 'a+') as fid_file:
		#~ fid_file.write('%.8g %.8g\n' % (karray[j], P1[j]))

plt.figure()
plt.plot(karray, Plin, label ='class')
plt.plot(karray, P1, label ='real space ps '+str(bias_model)+'-'+str(rsd))
plt.plot(karray, P2, label ='redshift space ps '+str(bias_model)+'-'+str(rsd))
plt.axvline(kmax, c='k', label='k max')
plt.title('The total neutrino mass is '+str(Mnu))
plt.legend(loc='upper right')
plt.xlabel('k [h/Mpc] ')
plt.ylabel(r'Pk $[h^3/ Mpc^{-3}]$')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()
