### code written by David Valcin
### for calibration details see arXiv:1901.06045

import camb
from camb import model, initialpower
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
kcase = '1'
# integration boundaries for k (in h/Mpc).
# It has to be within the boundaries of the k model above, respectively 0.15, 0.2, 0.4
kmin = 0.01
kmax = 0.15
kbins = 60
karray = np.logspace(np.log10(kmin), np.log10(kmax), kbins)
# desired redshift among the list [0.0, 0.5, 1.0, 2.0] for now. Interpolation might be added later
z = 0.5
# choice of neutrino mass between [0.0 - 0.6]eV
Mnu = 0.0
#choice of mass bins between (0,1,2,3) cf paper
mbin = 0
# choice of bias model between (1,2,3) cf paper: linear, polynomial, perturbation theory
bias_model = 1
#choice of RSD model between (1,2,3) cf paper
rsd = 1
# whether or not to include FoG (0 or 1), 1 for yes. sigma_v is defined with a default value fitted with the
# coefficients but the user can use different values
fog = 1
#sigma_v = [7.0, 7.0, 7.0, 7.0]



########################################################################
########################################################################
### Here an example 
### Import a linear power spectrum. Make sure to use values close the paper reference cosmology
########################################################################
########################################################################

#Now get matter power spectra and sigma8 at redshift 0 and 0.8
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(ns=0.9624, As = 2.13e-9)
#Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[z], kmax=2.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
kh, zk, Pk= results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)

#~ # interpolated the linear ps on desired k array
Plin = np.interp(karray, kh, Pk[0])


########################################################################
########################################################################
### Compute the non linear power spectrum
########################################################################
########################################################################



# A summary of your configuration is printed in the terminal
P1 = ps_calc(coord[0],kcase, Mnu, mbin, rsd, bias_model, karray, z, fog, Plin)
P2 = ps_calc(coord[1],kcase, Mnu, mbin, rsd, bias_model, karray, z, fog, Plin)


plt.figure()
plt.plot(karray, P1, label ='real space ps '+str(bias_model)+'-'+str(rsd))
plt.plot(karray, P2, label ='redshift space ps '+str(bias_model)+'-'+str(rsd))
plt.axvline(kmax, c='k', label='k max')
plt.title('The total neutrino mass is '+str(Mnu))
plt.legend(loc='upper right')
#~ plt.xlim(kmin, kmax)
#~ plt.ylim(1e2, 1e5)
plt.xlabel('k [h/Mpc] ')
plt.ylabel(r'Pk $[h^3/ Mpc^{-3}]$')
plt.xscale('log')
plt.yscale('log')
plt.show()
