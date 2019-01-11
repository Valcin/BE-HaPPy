import numpy as np
import h5py
import math
import readsnap
import matplotlib
#~ matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
sys.path.append('/home/david/codes/FAST-PT')
import myFASTPT as FPT
import scipy.interpolate as sp
import pyximport
pyximport.install()
import redshift_space_library as RSL
from readfof import FoF_catalog
import MAS_library as MASL
import Pk_library as PKL
import mass_function_library as MFL
import bias_library as BL
import tempfile
import expected_CF
import exp2
from table import table_write_pl
from time import time
from bias_library import halo_bias, bias
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.special import gamma
from fit_emcee import coeffit_pl,coeffit_pl2,coeffit_exp1, coeffit_exp2, coeffit_exp3,coeffit_Kaiser, coeffit_Scocci, coeffit_TNS, coeffit_eTNS


def poly(kstop, lb1, lb2, lb3, lb4, errlb1, errlb2, errlb3, errlb4, k, bias1,\
	bias2, bias3, bias4, errb1, errb2, errb3, errb4, mv, z, j, case):
		
	print ' the mass of neutrinos is '+str(mv)
	lim = np.where((k < kstop)&(k > 1e-2))[0]

	

	def funcb(k, b1, b2, b3, b4):
		return b1 + b2 * k**2 + b3 * k**3 + b4 * k**4 
		
	def funcbis(k, b1, b2, b4):
		return b1 + b2 * k**2 + b4 * k**4 
	
	
	#~ # here kh because the simu scale
	popF1, pcovF1 = curve_fit(funcb, k[lim], bias1[lim], sigma = errb1[lim],  check_finite=True, maxfev=500000)
	popF2, pcovF2 = curve_fit(funcb, k[lim], bias2[lim], sigma = errb2[lim],  check_finite=True, maxfev=500000)
	popF3, pcovF3 = curve_fit(funcb, k[lim], bias3[lim], sigma = errb3[lim],  check_finite=True, maxfev=500000)
	popF4, pcovF4 = curve_fit(funcb, k[lim], bias4[lim], sigma = errb4[lim],  check_finite=True, maxfev=500000)

	popF1bis, pcovF1bis = curve_fit(funcbis, k[lim], bias1[lim], sigma = errb1[lim],  check_finite=True, maxfev=500000)
	popF2bis, pcovF2bis = curve_fit(funcbis, k[lim], bias2[lim], sigma = errb2[lim],  check_finite=True, maxfev=500000)
	popF3bis, pcovF3bis = curve_fit(funcbis, k[lim], bias3[lim], sigma = errb3[lim],  check_finite=True, maxfev=500000)
	popF4bis, pcovF4bis = curve_fit(funcbis, k[lim], bias4[lim], sigma = errb4[lim],  check_finite=True, maxfev=500000)
		
		
		
	# odd power law----------------------------------------------------
	#~ b1x1_ml, b2x1_ml, b3x1_ml, b4x1_ml, b1x1_mcmc, b2x1_mcmc, b3x1_mcmc, b4x1_mcmc = coeffit_pl(kstop, lb1, errlb1, popF1, k, bias1, errb1)
	#~ b1x2_ml, b2x2_ml, b3x2_ml, b4x2_ml, b1x2_mcmc, b2x2_mcmc, b3x2_mcmc, b4x2_mcmc = coeffit_pl(kstop, lb2, errlb2, popF2, k, bias2, errb2)
	#~ b1x3_ml, b2x3_ml, b3x3_ml, b4x3_ml, b1x3_mcmc, b2x3_mcmc, b3x3_mcmc, b4x3_mcmc = coeffit_pl(kstop, lb3, errlb3, popF3, k, bias3, errb3)
	#~ b1x4_ml, b2x4_ml, b3x4_ml, b4x4_ml, b1x4_mcmc, b2x4_mcmc, b3x4_mcmc, b4x4_mcmc = coeffit_pl(kstop, lb4, errlb4, popF4, k, bias4, errb4)
	b1x1_ml, b2x1_ml, b3x1_ml, b4x1_ml = coeffit_pl(kstop, lb1, errlb1, popF1, k, bias1, errb1)
	b1x2_ml, b2x2_ml, b3x2_ml, b4x2_ml = coeffit_pl(kstop, lb2, errlb2, popF2, k, bias2, errb2)
	b1x3_ml, b2x3_ml, b3x3_ml, b4x3_ml = coeffit_pl(kstop, lb3, errlb3, popF3, k, bias3, errb3)
	b1x4_ml, b2x4_ml, b3x4_ml, b4x4_ml = coeffit_pl(kstop, lb4, errlb4, popF4, k, bias4, errb4)
	
	#-----
	#~ table_write_pl(b1x1_ml, b2x1_ml, b3x1_ml, b4x1_ml, b1x1_mcmc, b2x1_mcmc, b3x1_mcmc, b4x1_mcmc,\
	#~ b1x2_ml, b2x2_ml, b3x2_ml, b4x2_ml, b1x2_mcmc, b2x2_mcmc, b3x2_mcmc, b4x2_mcmc, \
	#~ b1x3_ml, b2x3_ml, b3x3_ml, b4x3_ml, b1x3_mcmc, b2x3_mcmc, b3x3_mcmc, b4x3_mcmc, \
	#~ b1x4_ml, b2x4_ml, b3x4_ml, b4x4_ml, b1x4_mcmc, b2x4_mcmc, b3x4_mcmc, b4x4_mcmc, 'coeffpl_'+str(z[j])+'.txt') 

	
	# even power law ----------------------------------------------------------------------------------------
	b1w1_ml, b2w1_ml, b4w1_ml = coeffit_pl2(kstop, lb1, errlb1, popF1bis, k, bias1, errb1)
	b1w2_ml, b2w2_ml, b4w2_ml = coeffit_pl2(kstop, lb2, errlb2, popF2bis, k, bias2, errb2)
	b1w3_ml, b2w3_ml, b4w3_ml = coeffit_pl2(kstop, lb3, errlb3, popF3bis, k, bias3, errb3)
	b1w4_ml, b2w4_ml, b4w4_ml = coeffit_pl2(kstop, lb4, errlb4, popF4bis, k, bias4, errb4)
		
		
#~ #### compute bias
#~ # power law odd ----------------------------------------------------------------------------
	#~ biasF1 = b1x1_mcmc[0] + b2x1_mcmc[0] * k**2 + b3x1_mcmc[0] * k**3 + b4x1_mcmc[0] * k**4
	#~ biasF2 = b1x2_mcmc[0] + b2x2_mcmc[0] * k**2 + b3x2_mcmc[0] * k**3 + b4x2_mcmc[0] * k**4
	#~ biasF3 = b1x3_mcmc[0] + b2x3_mcmc[0] * k**2 + b3x3_mcmc[0] * k**3 + b4x3_mcmc[0] * k**4
	#~ biasF4 = b1x4_mcmc[0] + b2x4_mcmc[0] * k**2 + b3x4_mcmc[0] * k**3 + b4x4_mcmc[0] * k**4
	biasF1 = b1x1_ml + b2x1_ml * k**2 + b3x1_ml * k**3 + b4x1_ml * k**4
	biasF2 = b1x2_ml + b2x2_ml * k**2 + b3x2_ml * k**3 + b4x2_ml * k**4
	biasF3 = b1x3_ml + b2x3_ml * k**2 + b3x3_ml * k**3 + b4x3_ml * k**4
	biasF4 = b1x4_ml + b2x4_ml * k**2 + b3x4_ml * k**3 + b4x4_ml * k**4

# power law even -------------------------------------------------------------------------------------------
	#~ biasF1bis = b1w1_mcmc[0] + b2w1_mcmc[0] * k**2 + b4w1_mcmc[0] * k**4
	#~ biasF2bis = b1w2_mcmc[0] + b2w2_mcmc[0] * k**2 + b4w2_mcmc[0] * k**4
	#~ biasF3bis = b1w3_mcmc[0] + b2w3_mcmc[0] * k**2 + b4w3_mcmc[0] * k**4
	#~ biasF4bis = b1w4_mcmc[0] + b2w4_mcmc[0] * k**2 + b4w4_mcmc[0] * k**4
	biasF1bis = b1w1_ml + b2w1_ml * k**2 + b4w1_ml * k**4
	biasF2bis = b1w2_ml + b2w2_ml * k**2 + b4w2_ml * k**4
	biasF3bis = b1w3_ml + b2w3_ml * k**2 + b4w3_ml * k**4
	biasF4bis = b1w4_ml + b2w4_ml * k**2 + b4w4_ml * k**4

#~ ###########################################################################
###########################################################################

	#~ cname1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/coeff_pl_'+str(mv)+'_z='+str(z[j])+'.txt'
	#~ cname1err = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/err_pl_'+str(mv)+'_z='+str(z[j])+'.txt'
	#~ cname1bis = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/coeff_ple_'+str(mv)+'_z='+str(z[j])+'.txt'
	#~ cname1errbis = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/err_ple_'+str(mv)+'_z='+str(z[j])+'.txt'
	
	cname1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/case'+str(case)+'/coeff_pl_'+str(mv)+'_z='+str(z[j])+'.txt'
	cname1err = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/case'+str(case)+'/err_pl_'+str(mv)+'_z='+str(z[j])+'.txt'
	cname1bis = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/case'+str(case)+'/coeff_ple_'+str(mv)+'_z='+str(z[j])+'.txt'
	cname1errbis = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(mv)+'eV/case'+str(case)+'/err_ple_'+str(mv)+'_z='+str(z[j])+'.txt'
		
	with open(cname1, 'w') as fid_file:
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (b1x1_ml, b2x1_ml, b3x1_ml, b4x1_ml))
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (b1x2_ml, b2x2_ml, b3x2_ml, b4x2_ml))
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (b1x3_ml, b2x3_ml, b3x3_ml, b4x3_ml))
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (b1x4_ml, b2x4_ml, b3x4_ml, b4x4_ml))
	fid_file.close()
	#~ with open(cname1err, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1x1_mcmc[1], b2x1_mcmc[1], b3x1_mcmc[1], b4x1_mcmc[1]\
		#~ ,b1x1_mcmc[2], b2x1_mcmc[2], b3x1_mcmc[2], b4x1_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1x2_mcmc[1], b2x2_mcmc[1], b3x2_mcmc[1], b4x2_mcmc[1]\
		#~ ,b1x2_mcmc[2], b2x2_mcmc[2], b3x2_mcmc[2], b4x2_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1x3_mcmc[1], b2x3_mcmc[1], b3x3_mcmc[1], b4x3_mcmc[1]\
		#~ ,b1x3_mcmc[2], b2x3_mcmc[2], b3x3_mcmc[2], b4x3_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1x4_mcmc[1], b2x4_mcmc[1], b3x4_mcmc[1], b4x4_mcmc[1]\
		#~ ,b1x4_mcmc[2], b2x4_mcmc[2], b3x4_mcmc[2], b4x4_mcmc[2]))
	#~ fid_file.close()
	with open(cname1bis, 'w') as fid_file:
		fid_file.write('%.8g %.8g %.8g\n' % (b1w1_ml, b2w1_ml, b4w1_ml))
		fid_file.write('%.8g %.8g %.8g\n' % (b1w2_ml, b2w2_ml, b4w2_ml))
		fid_file.write('%.8g %.8g %.8g\n' % (b1w3_ml, b2w3_ml, b4w3_ml))
		fid_file.write('%.8g %.8g %.8g\n' % (b1w4_ml, b2w4_ml, b4w4_ml))
	fid_file.close()
	#~ with open(cname1errbis, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1w1_mcmc[1], b2w1_mcmc[1], b4w1_mcmc[1]\
		#~ ,b1w1_mcmc[2], b2w1_mcmc[2], b4w1_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1w2_mcmc[1], b2w2_mcmc[1], b4w2_mcmc[1]\
		#~ ,b1w2_mcmc[2], b2w2_mcmc[2], b4w2_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1w3_mcmc[1], b2w3_mcmc[1], b4w3_mcmc[1]\
		#~ ,b1w3_mcmc[2], b2w3_mcmc[2], b4w3_mcmc[2]))
		#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g %.8g\n' % (b1w4_mcmc[1], b2w4_mcmc[1], b4w4_mcmc[1]\
		#~ ,b1w4_mcmc[2], b2w4_mcmc[2], b4w4_mcmc[2]))
	#~ fid_file.close()

##########################################################################

	#~ bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	#~ str(mv)+'eV/coeff_pl_'+str(mv)+'_z='+str(z[j])+'.txt')
	#~ bplbis = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	#~ str(mv)+'eV/coeff_ple_'+str(mv)+'_z='+str(z[j])+'.txt')
	
	#~ bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	#~ str(mv)+'eV/case'+str(case)+'/coeff_pl_'+str(mv)+'_z='+str(z[j])+'.txt')
	#~ bplbis =np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	#~ str(mv)+'eV/case'+str(case)+'/coeff_ple_'+str(mv)+'_z='+str(z[j])+'.txt')

	#~ b1pl = bpl[:,0]; b1plbis = bplbis[:,0]
	#~ b2pl = bpl[:,1]; b2plbis = bplbis[:,1]
	#~ b3pl = bpl[:,2]
	#~ b4pl = bpl[:,3]; b4plbis = bplbis[:,2]

#~ #### compute bias
#~ # power law odd ----------------------------------------------------------------------------
	#~ biasF1 = b1pl[0] + b2pl[0] * kbis**2 + b3pl[0] * kbis**3 + b4pl[0] * kbis**4
	#~ biasF2 = b1pl[1] + b2pl[1] * kbis**2 + b3pl[1] * kbis**3 + b4pl[1] * kbis**4
	#~ biasF3 = b1pl[2] + b2pl[2] * kbis**2 + b3pl[2] * kbis**3 + b4pl[2] * kbis**4
	#~ biasF4 = b1pl[3] + b2pl[3] * kbis**2 + b3pl[3] * kbis**3 + b4pl[3] * kbis**4

# power law even -------------------------------------------------------------------------------------------
	#~ biasF1bis = b1plbis[0] + b2plbis[0] * kbis**2 + b4plbis[0] * kbis**4
	#~ biasF2bis = b1plbis[1] + b2plbis[1] * kbis**2 + b4plbis[1] * kbis**4
	#~ biasF3bis = b1plbis[2] + b2plbis[2] * kbis**2 + b4plbis[2] * kbis**4
	#~ biasF4bis = b1plbis[3] + b2plbis[3] * kbis**2 + b4plbis[3] * kbis**4


##########################################################################

	return biasF1, biasF2, biasF3, biasF4, biasF1bis, biasF2bis, biasF3bis, biasF4bis
	#~ return biasF1bis, biasF2bis, biasF3bis, biasF4bis
	#~ return biasF1, biasF2, biasF3, biasF4
