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
from time import time
from load_data import ld_data
from loop_pt import pt_terms
from polynomial import poly
from perturbation import perturb
from bias_library import halo_bias, bias
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.special import gamma
from interp import interp_simu2
from fit_emcee import coeffit_pl,coeffit_pl2,coeffit_exp1, coeffit_exp2, coeffit_exp3,coeffit_Kaiser,coeffit_Kaiser2,\
coeffit_Scocci, coeffit_Scocci2, coeffit_TNS,coeffit_TNS2, coeffit_eTNS, coeffit_eTNS2

	
def RSD1(fz,fcc, Dz, j, kstop, kcamb, Pcamb, Pmm, biasF1, biasF2, biasF3, biasF4, k, Pmono1, Pmono2, Pmono3, \
		Pmono4, errPr1, errPr2, errPr3, errPr4, Pmod_dt, Pmod_tt, case,z,Mnu,A, B, C, D, E, F, G, H):
	####################################################################
	###### fit the Finger of God effect
	####################################################################
	#~ kbis = np.logspace(np.log10(np.min(kcamb)), np.log10(np.max(kcamb)), 200)
	#~ Plinbis = np.interp(kbis, k, Plin)

	#### compute tns coefficeints given mcmc results
	# set the parameters for the power spectrum window and
	# Fourier coefficient window 
	#P_window=np.array([.2,.2])  
	C_window=0.95

	# padding length 
	nu=-2; n_pad=len(kcamb)
	n_pad=int(0.5*len(kcamb))
	to_do=['all']
					
	# initialize the FASTPT class 
	# including extrapolation to higher and lower k  
	# time the operation
	t1=time()
	fastpt=FPT.FASTPT(kcamb,to_do=to_do,n_pad=n_pad, verbose=True) 
	t2=time()
	
####################################################################

	bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	str(Mnu)+'eV/case'+str(case)+'/coeff_pl_'+str(Mnu)+'_z='+str(z[j])+'.txt')
	b1pl = bpl[:,0]
	
	bpt3 = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	str(Mnu)+'eV/case'+str(case)+'/coeff_3exp_'+str(Mnu)+'_z='+str(z[j])+'.txt')
	b1pt3 = bpt3[:,0]
	b2pt3 = bpt3[:,1]
	bspt3 = bpt3[:,2]
	b3pt3 = bpt3[:,3]
	Npt3 = bpt3[:,4]
	
	AB2_1,AB4_1,AB6_1,AB8_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[0] ,C_window=C_window)
	AB2_2,AB4_2,AB6_2,AB8_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[1] ,C_window=C_window)
	AB2_3,AB4_3,AB6_3,AB8_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[2] ,C_window=C_window)
	AB2_4,AB4_4,AB6_4,AB8_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[3] ,C_window=C_window)
	
	
	AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[0] ,C_window=C_window)
	AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[1] ,C_window=C_window)
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[2] ,C_window=C_window)
	AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[3] ,C_window=C_window)

	AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
	AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = interp_simu2(z,j,k,kcamb, Pcamb,\
	AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
	AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4, 2)
		
		
	#~ dat_file_path = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/large_scale/'\
		#~ 'LS_z='+str(z[j])+'_.txt'
	#~ f = np.loadtxt(dat_file_path)
	#~ bls1 = f[0]
	#~ bls2 = f[1]
	#~ bls3 = f[2]
	#~ bls4 = f[3]
	#~ AB2ter_1,AB4ter_1,AB6ter_1,AB8ter_1 = fastpt.RSD_ABsum_components(Pcamb,fz,bls1 ,C_window=C_window)
	#~ AB2ter_2,AB4ter_2,AB6ter_2,AB8ter_2 = fastpt.RSD_ABsum_components(Pcamb,fz,bls2 ,C_window=C_window)
	#~ AB2ter_3,AB4ter_3,AB6ter_3,AB8ter_3 = fastpt.RSD_ABsum_components(Pcamb,fz,bls3 ,C_window=C_window)
	#~ AB2ter_4,AB4ter_4,AB6ter_4,AB8ter_4 = fastpt.RSD_ABsum_components(Pcamb,fz,bls4 ,C_window=C_window)
	
	#~ #-------------------------------------------------------
	#~ cname1m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl1_z='+str(z[j])+'.txt'
	#~ cname1m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl2_z='+str(z[j])+'.txt'
	#~ cname1m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl3_z='+str(z[j])+'.txt'
	#~ cname1m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl4_z='+str(z[j])+'.txt'
	#~ cname2m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt1_z='+str(z[j])+'.txt'
	#~ cname2m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt2_z='+str(z[j])+'.txt'
	#~ cname2m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt3_z='+str(z[j])+'.txt'
	#~ cname2m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt4_z='+str(z[j])+'.txt'
	#~ cname3m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin1_z='+str(z[j])+'.txt'
	#~ cname3m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin2_z='+str(z[j])+'.txt'
	#~ cname3m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin3_z='+str(z[j])+'.txt'
	#~ cname3m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin4_z='+str(z[j])+'.txt'


	#~ with open(cname1m1, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_1[index_k],AB4_1[index_k],AB6_1[index_k],AB8_1[index_k]))
	#~ with open(cname1m2, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_2[index_k],AB4_2[index_k],AB6_2[index_k],AB8_2[index_k]))
	#~ with open(cname1m3, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_3[index_k],AB4_3[index_k],AB6_3[index_k],AB8_3[index_k]))
	#~ with open(cname1m4, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_4[index_k],AB4_4[index_k],AB6_4[index_k],AB8_4[index_k]))
	#~ with open(cname2m1, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_1[index_k],AB4bis_1[index_k],AB6bis_1[index_k],AB8bis_1[index_k]))
	#~ with open(cname2m2, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_2[index_k],AB4bis_2[index_k],AB6bis_2[index_k],AB8bis_2[index_k]))
	#~ with open(cname2m3, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_3[index_k],AB4bis_3[index_k],AB6bis_3[index_k],AB8bis_3[index_k]))
	#~ with open(cname2m4, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_4[index_k],AB4bis_4[index_k],AB6bis_4[index_k],AB8bis_4[index_k]))
	#~ with open(cname3m1, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_1[index_k],AB4ter_1[index_k],AB6ter_1[index_k],AB8ter_1[index_k]))
	#~ with open(cname3m2, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_2[index_k],AB4ter_2[index_k],AB6ter_2[index_k],AB8ter_2[index_k]))
	#~ with open(cname3m3, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_3[index_k],AB4ter_3[index_k],AB6ter_3[index_k],AB8ter_3[index_k]))
	#~ with open(cname3m4, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_4[index_k],AB4ter_4[index_k],AB6ter_4[index_k],AB8ter_4[index_k]))

	

	#----------------------------------------------------------------
	#### compute mcmc coefficient of halo ps fit
	print 'kaiser'
	#~ b1,b2,b3,b4,bK1 = coeffit_Kaiser2(j, fcc, kstop,Pmm, biasF1, k, Pmono1, errPr1)
	bK1 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF1, k, Pmono1, errPr1)
	bK2 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF2, k, Pmono2, errPr2)
	bK3 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF3, k, Pmono3, errPr3)
	bK4 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF4, k, Pmono4, errPr4)
	#~ cn1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt'
	#~ with open(cn1, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bK1[0][0],bK2[0][0],bK3[0][0],bK4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Scoccimaro'
	bsco1 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF1, k, Pmono1, errPr1)
	bsco2 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF2, k, Pmono2, errPr2)
	bsco3 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF3, k, Pmono3, errPr3)
	bsco4 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF4, k, Pmono4, errPr4)
	#~ cn2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt'
	#~ with open(cn2, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bsco1[0][0],bsco2[0][0],bsco3[0][0],bsco4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Tns'
	btns1 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF1, k, Pmono1, errPr1, AB2_1, AB4_1, AB6_1, AB8_1)
	btns2 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF2, k, Pmono2, errPr2, AB2_2, AB4_2, AB6_2, AB8_2)
	btns3 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF3, k, Pmono3, errPr3, AB2_3, AB4_3, AB6_3, AB8_3)
	btns4 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF4, k, Pmono4, errPr4, AB2_4, AB4_4, AB6_4, AB8_4)
	#~ cn3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt'
	#~ with open(cn3, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (btns1[0][0],btns2[0][0],btns3[0][0],btns4[0][0]))
	#~ fid_file.close()

	#~ #----------------------------------------------------------------------------------------
	print 'eTns'
	betns1 = coeffit_eTNS(j, fcc, kstop, b1pt3[0], b2pt3[0], bspt3[0], b3pt3[0], Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono1, errPr1, AB2bis_1, AB4bis_1,\
	AB6bis_1, AB8bis_1, Npt3[0])
	betns2 = coeffit_eTNS(j, fcc, kstop, b1pt3[1], b2pt3[1], bspt3[1], b3pt3[1], Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono2, errPr2, AB2bis_2, AB4bis_2,\
	AB6bis_2, AB8bis_2, Npt3[1])
	betns3 = coeffit_eTNS(j, fcc, kstop, b1pt3[2], b2pt3[2], bspt3[2], b3pt3[2], Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono3, errPr3, AB2bis_3, AB4bis_3,\
	AB6bis_3, AB8bis_3, Npt3[2])
	betns4 = coeffit_eTNS(j, fcc, kstop, b1pt3[3], b2pt3[3], bspt3[3], b3pt3[3], Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono4, errPr4, AB2bis_4, AB4bis_4,\
	AB6bis_4, AB8bis_4, Npt3[3])
	#~ cn4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt'
	#~ with open(cn4, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (betns1[0][0],betns2[0][0],betns3[0][0],betns4[0][0]))
	#~ fid_file.close()
	
##########################################################################
##########################################################################

	#### compute the different power spectra given the mcmc results
	bK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')
	#~ BK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkaibis_z='+str(z[j])+'.txt')
	#~ biasK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')

	def kaips(b,sigma):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return Pmm*(b**2*coeffA +  2/3.*b*fcc*coeffB + 1/5.*fcc**2*coeffC)
		#~ return Pmm*(b**2 +  2/3.*b*fcc + 1/5.*fcc**2)
		
	#~ kai1 = kaips(biasF1, bK[0])
	#~ kai2 = kaips(biasF2, bK[1])
	#~ kai3 = kaips(biasF3, bK[2])
	#~ kai4 = kaips(biasF4, bK[3])
	kai1 = kaips(biasF1, bK1)
	kai2 = kaips(biasF2, bK2)
	kai3 = kaips(biasF3, bK3)
	kai4 = kaips(biasF4, bK4)
	#*****************************
	#~ k1 = kaips(bF1, BK[0])
	#~ k2 = kaips(bF2, BK[1])
	#~ k3 = kaips(bF3, BK[2])
	#~ k4 = kaips(bF4, BK[3])

	#~ k1ter = kaips(bF1, biasK[0])
	#~ k2ter = kaips(bF2, biasK[1])
	#~ k3ter = kaips(bF3, biasK[2])
	#~ k4ter = kaips(bF4, biasK[3])

	#---------------------------------------------------------------------------------------
	bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	#~ Bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispscobis_z='+str(z[j])+'.txt')
	#~ bs = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	def scops(b,sigma):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return b**2*Pmm*coeffA + 2/3.*b*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt

	#~ sco1 = scops(biasF1, bsco[0])
	#~ sco2 = scops(biasF2, bsco[1])
	#~ sco3 = scops(biasF3, bsco[2])
	#~ sco4 = scops(biasF4, bsco[3])
	sco1 = scops(biasF1, bsco1)
	sco2 = scops(biasF2, bsco2)
	sco3 = scops(biasF3, bsco3)
	sco4 = scops(biasF4, bsco4)
	
	#~ s1 = scops(bF1, Bsco[0])
	#~ s2 = scops(bF2, Bsco[1])
	#~ s3 = scops(bF3, Bsco[2])
	#~ s4 = scops(bF4, Bsco[3])
	
	#~ s1ter = scops(bF1, bs[0])
	#~ s2ter = scops(bF2, bs[1])
	#~ s3ter = scops(bF3, bs[2])
	#~ s4ter = scops(bF4, bs[3])
	#~ #---------------------------------------------------------------------------------------
	btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	#~ Btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptnsbis_z='+str(z[j])+'.txt')
	#~ bt = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	def tnsps(b,sigma, AB2, AB4, AB6, AB8):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		return b**2*Pmm*coeffA + 2/3.*b*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
		+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8)

	#~ tns1 = tnsps(biasF1,btns[0], AB2_1, AB4_1, AB6_1, AB8_1)
	#~ tns2 = tnsps(biasF2,btns[1], AB2_2, AB4_2, AB6_2, AB8_2)
	#~ tns3 = tnsps(biasF3,btns[2], AB2_3, AB4_3, AB6_3, AB8_3)
	#~ tns4 = tnsps(biasF4,btns[3], AB2_4, AB4_4, AB6_4, AB8_4)
	tns1 = tnsps(biasF1,btns1, AB2_1, AB4_1, AB6_1, AB8_1)
	tns2 = tnsps(biasF2,btns2, AB2_2, AB4_2, AB6_2, AB8_2)
	tns3 = tnsps(biasF3,btns3, AB2_3, AB4_3, AB6_3, AB8_3)
	tns4 = tnsps(biasF4,btns4, AB2_4, AB4_4, AB6_4, AB8_4)
	
	#~ t1 = tnsps(bF1,Btns[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2 = tnsps(bF2,Btns[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3 = tnsps(bF3,Btns[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4 = tnsps(bF4,Btns[3], ab2_4, ab4_4, ab6_4, ab8_4)
	
	#~ t1ter = tnsps(bF1,bt[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2ter = tnsps(bF2,bt[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3ter = tnsps(bF3,bt[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4ter = tnsps(bF4,bt[3], ab2_4, ab4_4, ab6_4, ab8_4)
	#-------------------------------------------------------------------
	betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	#~ Betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetnsbis_z='+str(z[j])+'.txt')
	#~ be = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	def etnsps(b1,b2,bs,b3nl,sigma, AB2, AB4, AB6, AB8, N, sca=None):
		PsptD1z = b1**2*Pmm + b1*b2*A+ 1/4.*b2**2*B+ b1*bs*C+ 1/2.*b2*bs*D+ 1/4.*bs**2*E+ 2*b1*b3nl*F + N
		PsptT = b1* Pmod_dt+ b2*G+ bs*H + b3nl*F
		#~ kappa = k*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		if sca:
			return  PsptD1z*coeffA*sca**2 + 2/3.*fcc*PsptT*coeffB*sca + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		else:
			return  PsptD1z*coeffA + 2/3.*fcc*PsptT*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return  PsptD1z + 2/3.*fcc*PsptT + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
		
	
	#~ etns1 = etnsps(b1pt3[0], b2pt3[0], bspt3[0], b3pt3[0], betns[0], AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1)  
	#~ etns2 = etnsps(b1pt3[1], b2pt3[1], bspt3[1], b3pt3[1], betns[1], AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2)  
	#~ etns3 = etnsps(b1pt3[2], b2pt3[2], bspt3[2], b3pt3[2], betns[2], AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3)  
	#~ etns4 = etnsps(b1pt3[3], b2pt3[3], bspt3[3], b3pt3[3], betns[3], AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4) 
	etns1 = etnsps(b1pt3[0], b2pt3[0], bspt3[0], b3pt3[0], betns1, AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1, Npt3[0])  
	etns2 = etnsps(b1pt3[1], b2pt3[1], bspt3[1], b3pt3[1], betns2, AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2, Npt3[1])  
	etns3 = etnsps(b1pt3[2], b2pt3[2], bspt3[2], b3pt3[2], betns3, AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3, Npt3[2])  
	etns4 = etnsps(b1pt3[3], b2pt3[3], bspt3[3], b3pt3[3], betns4, AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4, Npt3[3]) 
	 
	#~ e1 = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], Betns[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1,(Bias_eff_t1/Bias_eff0_t1))  
	#~ e2 = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], Betns[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2,(Bias_eff_t2/Bias_eff0_t2))  
	#~ e3 = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], Betns[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3,(Bias_eff_t3/Bias_eff0_t3))  
	#~ e4 = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], Betns[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4,(Bias_eff_t4/Bias_eff0_t4))
	  
	#~ e1ter = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], be[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1)  
	#~ e2ter = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], be[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2)  
	#~ e3ter = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], be[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3)  
	#~ e4ter = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], be[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4)  

	#~ return kai1, kai2, kai3, kai4, 
	#~ return sco1, sco2, sco3, sco4
	#~ return tns1, tns2, tns3, tns4
	#~ return etns1, etns2, etns3, etns4
	return kai1, kai2, kai3, kai4, sco1, sco2, sco3, sco4, tns1, tns2, tns3, tns4, etns1, etns2, etns3, etns4
	
################################################################################################################
################################################################################################################
################################################################################################################
	
def RSD2(fz,fcc, Dz, j, kstop, kcamb, Pcamb, Pmm, biasF1, biasF2, biasF3, biasF4, k, Pmono1, Pmono2, Pmono3, \
		Pmono4, errPr1, errPr2, errPr3, errPr4, Pmod_dt, Pmod_tt, case,z,Mnu,A, B, C, D, E, F, G, H, sca1,
		sca2, sca3, sca4 ):
	####################################################################
	###### fit the Finger of God effect
	####################################################################
	#~ kbis = np.logspace(np.log10(np.min(kcamb)), np.log10(np.max(kcamb)), 200)
	#~ Plinbis = np.interp(kbis, k, Plin)

	#### compute tns coefficeints given mcmc results
	# set the parameters for the power spectrum window and
	# Fourier coefficient window 
	#P_window=np.array([.2,.2])  
	C_window=0.95

	# padding length 
	nu=-2; n_pad=len(kcamb)
	n_pad=int(0.5*len(kcamb))
	to_do=['all']
					
	# initialize the FASTPT class 
	# including extrapolation to higher and lower k  
	# time the operation
	t1=time()
	fastpt=FPT.FASTPT(kcamb,to_do=to_do,n_pad=n_pad, verbose=True) 
	t2=time()
	
####################################################################
	if sca1 is not None:
		Mnu0 = 0.0
		bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
		str(Mnu0)+'eV/case'+str(case)+'/coeff_pl_'+str(Mnu0)+'_z='+str(z[j])+'.txt')
		b1pl = bpl[:,0]
		
		bpt3 = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
		str(Mnu0)+'eV/case'+str(case)+'/coeff_3exp_'+str(Mnu0)+'_z='+str(z[j])+'.txt')
		b1pt3 = bpt3[:,0]
		b2pt3 = bpt3[:,1]
		bspt3 = bpt3[:,2]
		b3pt3 = bpt3[:,3]
		Npt3 = bpt3[:,4]
		
		AB2_1,AB4_1,AB6_1,AB8_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[0]*sca1 ,C_window=C_window)
		AB2_2,AB4_2,AB6_2,AB8_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[1]*sca2 ,C_window=C_window)
		AB2_3,AB4_3,AB6_3,AB8_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[2]*sca3 ,C_window=C_window)
		AB2_4,AB4_4,AB6_4,AB8_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[3]*sca4 ,C_window=C_window)
		
		AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[0]*sca1 ,C_window=C_window)
		AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[1]*sca2 ,C_window=C_window)
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[2]*sca3 ,C_window=C_window)
		AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[3]*sca4 ,C_window=C_window)
		
		AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
		AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = interp_simu2(z,j, k,kcamb, Pcamb,\
		AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
		AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4, 2)
		
	else:
		bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
		str(Mnu)+'eV/case'+str(case)+'/coeff_pl_'+str(Mnu)+'_z='+str(z[j])+'.txt')
		b1pl = bpl[:,0]
		
		bpt3 = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
		str(Mnu)+'eV/case'+str(case)+'/coeff_3exp_'+str(Mnu)+'_z='+str(z[j])+'.txt')
		b1pt3 = bpt3[:,0]
		b2pt3 = bpt3[:,1]
		bspt3 = bpt3[:,2]
		b3pt3 = bpt3[:,3]
		Npt3 = bpt3[:,4]
		
		AB2_1,AB4_1,AB6_1,AB8_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[0] ,C_window=C_window)
		AB2_2,AB4_2,AB6_2,AB8_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[1] ,C_window=C_window)
		AB2_3,AB4_3,AB6_3,AB8_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[2] ,C_window=C_window)
		AB2_4,AB4_4,AB6_4,AB8_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[3] ,C_window=C_window)
		
		
		AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[0] ,C_window=C_window)
		AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[1] ,C_window=C_window)
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[2] ,C_window=C_window)
		AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[3] ,C_window=C_window)
	
		AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
		AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = interp_simu2(z,j,k,kcamb, Pcamb,\
		AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
		AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
		AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4, 2)
		
		
		
		
	#~ dat_file_path = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/large_scale/'\
		#~ 'LS_z='+str(z[j])+'_.txt'
	#~ f = np.loadtxt(dat_file_path)
	#~ bls1 = f[0]
	#~ bls2 = f[1]
	#~ bls3 = f[2]
	#~ bls4 = f[3]
	#~ AB2ter_1,AB4ter_1,AB6ter_1,AB8ter_1 = fastpt.RSD_ABsum_components(Pcamb,fz,bls1 ,C_window=C_window)
	#~ AB2ter_2,AB4ter_2,AB6ter_2,AB8ter_2 = fastpt.RSD_ABsum_components(Pcamb,fz,bls2 ,C_window=C_window)
	#~ AB2ter_3,AB4ter_3,AB6ter_3,AB8ter_3 = fastpt.RSD_ABsum_components(Pcamb,fz,bls3 ,C_window=C_window)
	#~ AB2ter_4,AB4ter_4,AB6ter_4,AB8ter_4 = fastpt.RSD_ABsum_components(Pcamb,fz,bls4 ,C_window=C_window)
	
	#~ #-------------------------------------------------------
	#~ cname1m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl1_z='+str(z[j])+'.txt'
	#~ cname1m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl2_z='+str(z[j])+'.txt'
	#~ cname1m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl3_z='+str(z[j])+'.txt'
	#~ cname1m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pl4_z='+str(z[j])+'.txt'
	#~ cname2m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt1_z='+str(z[j])+'.txt'
	#~ cname2m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt2_z='+str(z[j])+'.txt'
	#~ cname2m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt3_z='+str(z[j])+'.txt'
	#~ cname2m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_pt4_z='+str(z[j])+'.txt'
	#~ cname3m1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin1_z='+str(z[j])+'.txt'
	#~ cname3m2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin2_z='+str(z[j])+'.txt'
	#~ cname3m3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin3_z='+str(z[j])+'.txt'
	#~ cname3m4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/TNS_coeff/'\
	#~ 'tns_lin4_z='+str(z[j])+'.txt'


	#~ with open(cname1m1, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_1[index_k],AB4_1[index_k],AB6_1[index_k],AB8_1[index_k]))
	#~ with open(cname1m2, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_2[index_k],AB4_2[index_k],AB6_2[index_k],AB8_2[index_k]))
	#~ with open(cname1m3, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_3[index_k],AB4_3[index_k],AB6_3[index_k],AB8_3[index_k]))
	#~ with open(cname1m4, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2_4[index_k],AB4_4[index_k],AB6_4[index_k],AB8_4[index_k]))
	#~ with open(cname2m1, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_1[index_k],AB4bis_1[index_k],AB6bis_1[index_k],AB8bis_1[index_k]))
	#~ with open(cname2m2, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_2[index_k],AB4bis_2[index_k],AB6bis_2[index_k],AB8bis_2[index_k]))
	#~ with open(cname2m3, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_3[index_k],AB4bis_3[index_k],AB6bis_3[index_k],AB8bis_3[index_k]))
	#~ with open(cname2m4, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2bis_4[index_k],AB4bis_4[index_k],AB6bis_4[index_k],AB8bis_4[index_k]))
	#~ with open(cname3m1, 'w') as fid_file:
		#~ for index_k in xrange(len(kbis)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_1[index_k],AB4ter_1[index_k],AB6ter_1[index_k],AB8ter_1[index_k]))
	#~ with open(cname3m2, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_2[index_k],AB4ter_2[index_k],AB6ter_2[index_k],AB8ter_2[index_k]))
	#~ with open(cname3m3, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_3[index_k],AB4ter_3[index_k],AB6ter_3[index_k],AB8ter_3[index_k]))
	#~ with open(cname3m4, 'w') as fid_file:
		#~ for index_k in xrange(len(k)):
			#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (AB2ter_4[index_k],AB4ter_4[index_k],AB6ter_4[index_k],AB8ter_4[index_k]))

	

	#----------------------------------------------------------------
	#### compute mcmc coefficient of halo ps fit
	print 'kaiser'
	bk1a,bk2a,bk3a,bk4a,bK1 = coeffit_Kaiser2(j, fcc, kstop,Pmm, k, Pmono1, errPr1)
	bk1b,bk2b,bk3b,bk4b,bK2 = coeffit_Kaiser2(j, fcc, kstop,Pmm, k, Pmono2, errPr2)
	bk1c,bk2c,bk3c,bk4c,bK3 = coeffit_Kaiser2(j, fcc, kstop,Pmm, k, Pmono3, errPr3)
	bk1d,bk2d,bk3d,bk4d,bK4 = coeffit_Kaiser2(j, fcc, kstop,Pmm, k, Pmono4, errPr4)
	
	#~ cn1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt'
	#~ with open(cn1, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bK1[0][0],bK2[0][0],bK3[0][0],bK4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Scoccimaro'
	bs1a,bs2a,bs3a,bs4a,bsco1 = coeffit_Scocci2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono1, errPr1)
	bs1b,bs2b,bs3b,bs4b,bsco2 = coeffit_Scocci2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono2, errPr2)
	bs1c,bs2c,bs3c,bs4c,bsco3 = coeffit_Scocci2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono3, errPr3)
	bs1d,bs2d,bs3d,bs4d,bsco4 = coeffit_Scocci2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono4, errPr4)
	#~ cn2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt'
	#~ with open(cn2, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bsco1[0][0],bsco2[0][0],bsco3[0][0],bsco4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Tns'
	bt1a,bt2a,bt3a,bt4a, btns1 = coeffit_TNS2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono1, errPr1, AB2_1, AB4_1, AB6_1, AB8_1)
	bt1b,bt2b,bt3b,bt4b, btns2 = coeffit_TNS2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono2, errPr2, AB2_2, AB4_2, AB6_2, AB8_2)
	bt1c,bt2c,bt3c,bt4c, btns3 = coeffit_TNS2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono3, errPr3, AB2_3, AB4_3, AB6_3, AB8_3)
	bt1d,bt2d,bt3d,bt4d, btns4 = coeffit_TNS2(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, k, Pmono4, errPr4, AB2_4, AB4_4, AB6_4, AB8_4)
	#~ cn3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt'
	#~ with open(cn3, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (btns1[0][0],btns2[0][0],btns3[0][0],btns4[0][0]))
	#~ fid_file.close()

	#~ #----------------------------------------------------------------------------------------
	print 'eTns'
	if sca1 is not None:
		be1a, be2a, besa, be3a, Nea, betns1 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono1, errPr1, AB2bis_1, AB4bis_1,\
		AB6bis_1, AB8bis_1, Npt3[0])
		be1b, be2b, besb, be3b, Neb, betns2 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono2, errPr2, AB2bis_2, AB4bis_2,\
		AB6bis_2, AB8bis_2, Npt3[1])
		be1c, be2c, besc, be3c, Nec, betns3 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono3, errPr3, AB2bis_3, AB4bis_3,\
		AB6bis_3, AB8bis_3, Npt3[2])
		be1d, be2d, besd, be3d, Ned, betns4 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono4, errPr4, AB2bis_4, AB4bis_4,\
		AB6bis_4, AB8bis_4, Npt3[3])
	else:
		be1a, be2a, Nea, betns1 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono1, errPr1, AB2bis_1, AB4bis_1,\
		AB6bis_1, AB8bis_1, Npt3[0])
		be1b, be2b, Neb, betns2 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono2, errPr2, AB2bis_2, AB4bis_2,\
		AB6bis_2, AB8bis_2, Npt3[1])
		be1c, be2c, Nec, betns3 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono3, errPr3, AB2bis_3, AB4bis_3,\
		AB6bis_3, AB8bis_3, Npt3[2])
		be1d, be2d, Ned, betns4 = coeffit_eTNS2(j, fcc, kstop, Pmm, Pmod_dt, Pmod_tt,\
		A, B, C, D, E, F, G, H, k, Pmono4, errPr4, AB2bis_4, AB4bis_4,\
		AB6bis_4, AB8bis_4, Npt3[3])
	#~ cn4 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt'
	#~ with open(cn4, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (betns1[0][0],betns2[0][0],betns3[0][0],betns4[0][0]))
	#~ fid_file.close()
	
##########################################################################
##########################################################################

	#### compute the different power spectra given the mcmc results
	bK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')
	#~ BK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkaibis_z='+str(z[j])+'.txt')
	#~ biasK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')

	def kaips(b1,b2,b3,b4,sigma):
		bcc = b1 + b2*k**2 + b3*k**3 + b4*k**4 
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return Pmm*(bcc**2*coeffA +  2/3.*bcc*fcc*coeffB + 1/5.*fcc**2*coeffC)
		#~ return Pmm*(b**2 +  2/3.*b*fcc + 1/5.*fcc**2)
		
	#~ kai1 = kaips(biasF1, bK[0])
	#~ kai2 = kaips(biasF2, bK[1])
	#~ kai3 = kaips(biasF3, bK[2])
	#~ kai4 = kaips(biasF4, bK[3])
	kai1 = kaips(bk1a,bk2a,bk3a,bk4a, bK1)
	kai2 = kaips(bk1b,bk2b,bk3b,bk4b, bK2)
	kai3 = kaips(bk1c,bk2c,bk3c,bk4c, bK3)
	kai4 = kaips(bk1d,bk2d,bk3d,bk4d, bK4)
	#*****************************
	#~ k1 = kaips(bF1, BK[0])
	#~ k2 = kaips(bF2, BK[1])
	#~ k3 = kaips(bF3, BK[2])
	#~ k4 = kaips(bF4, BK[3])

	#~ k1ter = kaips(bF1, biasK[0])
	#~ k2ter = kaips(bF2, biasK[1])
	#~ k3ter = kaips(bF3, biasK[2])
	#~ k4ter = kaips(bF4, biasK[3])

	#---------------------------------------------------------------------------------------
	bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	#~ Bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispscobis_z='+str(z[j])+'.txt')
	#~ bs = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	def scops(b1,b2,b3,b4,sigma):
		bcc = b1 + b2*k**2 + b3*k**3 + b4*k**4 
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return bcc**2*Pmm*coeffA + 2/3.*bcc*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt

	#~ sco1 = scops(biasF1, bsco[0])
	#~ sco2 = scops(biasF2, bsco[1])
	#~ sco3 = scops(biasF3, bsco[2])
	#~ sco4 = scops(biasF4, bsco[3])
	sco1 = scops(bs1a,bs2a,bs3a,bs4a, bsco1)
	sco2 = scops(bs1b,bs2b,bs3b,bs4b, bsco2)
	sco3 = scops(bs1c,bs2c,bs3c,bs4c, bsco3)
	sco4 = scops(bs1d,bs2d,bs3d,bs4d, bsco4)
	
	#~ s1 = scops(bF1, Bsco[0])
	#~ s2 = scops(bF2, Bsco[1])
	#~ s3 = scops(bF3, Bsco[2])
	#~ s4 = scops(bF4, Bsco[3])
	
	#~ s1ter = scops(bF1, bs[0])
	#~ s2ter = scops(bF2, bs[1])
	#~ s3ter = scops(bF3, bs[2])
	#~ s4ter = scops(bF4, bs[3])
	#~ #---------------------------------------------------------------------------------------
	btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	#~ Btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptnsbis_z='+str(z[j])+'.txt')
	#~ bt = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	def tnsps(b1,b2,b3,b4,sigma, AB2, AB4, AB6, AB8):
		bcc =  b1 + b2*k**2 + b3*k**3 + b4*k**4 
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		return bcc**2*Pmm*coeffA + 2/3.*bcc*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
		+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8)

	#~ tns1 = tnsps(biasF1,btns[0], AB2_1, AB4_1, AB6_1, AB8_1)
	#~ tns2 = tnsps(biasF2,btns[1], AB2_2, AB4_2, AB6_2, AB8_2)
	#~ tns3 = tnsps(biasF3,btns[2], AB2_3, AB4_3, AB6_3, AB8_3)
	#~ tns4 = tnsps(biasF4,btns[3], AB2_4, AB4_4, AB6_4, AB8_4)
	tns1 = tnsps(bt1a,bt2a,bt3a,bt4a,btns1, AB2_1, AB4_1, AB6_1, AB8_1)
	tns2 = tnsps(bt1b,bt2b,bt3b,bt4b,btns2, AB2_2, AB4_2, AB6_2, AB8_2)
	tns3 = tnsps(bt1c,bt2c,bt3c,bt4c,btns3, AB2_3, AB4_3, AB6_3, AB8_3)
	tns4 = tnsps(bt1d,bt2d,bt3d,bt4d,btns4, AB2_4, AB4_4, AB6_4, AB8_4)
	
	#~ t1 = tnsps(bF1,Btns[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2 = tnsps(bF2,Btns[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3 = tnsps(bF3,Btns[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4 = tnsps(bF4,Btns[3], ab2_4, ab4_4, ab6_4, ab8_4)
	
	#~ t1ter = tnsps(bF1,bt[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2ter = tnsps(bF2,bt[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3ter = tnsps(bF3,bt[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4ter = tnsps(bF4,bt[3], ab2_4, ab4_4, ab6_4, ab8_4)
	#-------------------------------------------------------------------
	betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	#~ Betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetnsbis_z='+str(z[j])+'.txt')
	#~ be = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	def etnsps(b1,b2,sigma, AB2, AB4, AB6, AB8, N, sca=None):
		bs = -4/7.*(b1-1)
		b3nl = 32/315.*(b1-1)
		PsptD1z = b1**2*Pmm + b1*b2*A+ 1/4.*b2**2*B+ b1*bs*C+ 1/2.*b2*bs*D+ 1/4.*bs**2*E+ 2*b1*b3nl*F + N
		PsptT = b1* Pmod_dt+ b2*G+ bs*H + b3nl*F
		#~ kappa = k*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		if sca:
			return  PsptD1z*coeffA*sca**2 + 2/3.*fcc*PsptT*coeffB*sca + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		else:
			return  PsptD1z*coeffA + 2/3.*fcc*PsptT*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return  PsptD1z + 2/3.*fcc*PsptT + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
		
	
	#~ etns1 = etnsps(b1pt3[0], b2pt3[0], bspt3[0], b3pt3[0], betns[0], AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1)  
	#~ etns2 = etnsps(b1pt3[1], b2pt3[1], bspt3[1], b3pt3[1], betns[1], AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2)  
	#~ etns3 = etnsps(b1pt3[2], b2pt3[2], bspt3[2], b3pt3[2], betns[2], AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3)  
	#~ etns4 = etnsps(b1pt3[3], b2pt3[3], bspt3[3], b3pt3[3], betns[3], AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4) 
	etns1 = etnsps(be1a, be2a,  betns1, AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1, Nea)  
	etns2 = etnsps(be1b, be2b,  betns2, AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2, Neb)  
	etns3 = etnsps(be1c, be2c,  betns3, AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3, Nec)  
	etns4 = etnsps(be1d, be2d,  betns4, AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4, Ned) 
	 
	#~ e1 = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], Betns[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1,(Bias_eff_t1/Bias_eff0_t1))  
	#~ e2 = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], Betns[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2,(Bias_eff_t2/Bias_eff0_t2))  
	#~ e3 = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], Betns[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3,(Bias_eff_t3/Bias_eff0_t3))  
	#~ e4 = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], Betns[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4,(Bias_eff_t4/Bias_eff0_t4))
	  
	#~ e1ter = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], be[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1)  
	#~ e2ter = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], be[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2)  
	#~ e3ter = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], be[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3)  
	#~ e4ter = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], be[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4)  

	#~ return kai1, kai2, kai3, kai4, 
	#~ return sco1, sco2, sco3, sco4
	#~ return tns1, tns2, tns3, tns4
	#~ return etns1, etns2, etns3, etns4
	return kai1, kai2, kai3, kai4, sco1, sco2, sco3, sco4, tns1, tns2, tns3, tns4, etns1, etns2, etns3, etns4
	
	
def RSD3(fz,fcc, Dz, j, kstop, kcamb, Pcamb, Pmm, biasF1, biasF2, biasF3, biasF4, k, Pmono1, Pmono2, Pmono3, \
		Pmono4, errPr1, errPr2, errPr3, errPr4, Pmod_dt, Pmod_tt, case,z,Mnu,A, B, C, D, E, F, G, H, sca1 = None,
		sca2 = None, sca3 = None, sca4 = None ):
	####################################################################
	###### fit the Finger of God effect
	####################################################################
	#~ kbis = np.logspace(np.log10(np.min(kcamb)), np.log10(np.max(kcamb)), 200)
	#~ Plinbis = np.interp(kbis, k, Plin)

	#### compute tns coefficeints given mcmc results
	# set the parameters for the power spectrum window and
	# Fourier coefficient window 
	#P_window=np.array([.2,.2])  
	C_window=0.95

	# padding length 
	nu=-2; n_pad=len(kcamb)
	n_pad=int(0.5*len(kcamb))
	to_do=['all']
					
	# initialize the FASTPT class 
	# including extrapolation to higher and lower k  
	# time the operation
	t1=time()
	fastpt=FPT.FASTPT(kcamb,to_do=to_do,n_pad=n_pad, verbose=True) 
	t2=time()
	
####################################################################

	Mnu0 = 0.0
	bpl = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	str(Mnu0)+'eV/case'+str(case)+'/coeff_pl_'+str(Mnu0)+'_z='+str(z[j])+'.txt')
	b1pl = bpl[:,0]
	
	bpt3 = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+\
	str(Mnu0)+'eV/case'+str(case)+'/coeff_3exp_'+str(Mnu0)+'_z='+str(z[j])+'.txt')
	b1pt3 = bpt3[:,0]
	b2pt3 = bpt3[:,1]
	bspt3 = bpt3[:,2]
	b3pt3 = bpt3[:,3]
	Npt3 = bpt3[:,4]
	
	AB2_1,AB4_1,AB6_1,AB8_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[0]*sca1 ,C_window=C_window)
	AB2_2,AB4_2,AB6_2,AB8_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[1]*sca2 ,C_window=C_window)
	AB2_3,AB4_3,AB6_3,AB8_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[2]*sca3 ,C_window=C_window)
	AB2_4,AB4_4,AB6_4,AB8_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pl[3]*sca4 ,C_window=C_window)
	
	AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[0]*sca1 ,C_window=C_window)
	AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[1]*sca2 ,C_window=C_window)
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[2]*sca3 ,C_window=C_window)
	AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = fastpt.RSD_ABsum_components(Pcamb,fz,b1pt3[3]*sca4 ,C_window=C_window)
	
	AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
	AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4 = interp_simu2(z,j, k,kcamb, Pcamb,\
	AB2_1,AB4_1,AB6_1,AB8_1,AB2_2,AB4_2,AB6_2,AB8_2,AB2_3,AB4_3,AB6_3,AB8_3,\
	AB2_4,AB4_4,AB6_4,AB8_4,AB2bis_1,AB4bis_1,AB6bis_1,AB8bis_1 ,AB2bis_2,AB4bis_2,AB6bis_2,AB8bis_2,\
	AB2bis_3,AB4bis_3,AB6bis_3,AB8bis_3, AB2bis_4,AB4bis_4,AB6bis_4,AB8bis_4, 2)
		


	#----------------------------------------------------------------
	#### compute mcmc coefficient of halo ps fit
	print 'kaiser'
	#~ b1,b2,b3,b4,bK1 = coeffit_Kaiser2(j, fcc, kstop,Pmm, biasF1, k, Pmono1, errPr1)
	bK1 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF1, k, Pmono1, errPr1)
	bK2 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF2, k, Pmono2, errPr2)
	bK3 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF3, k, Pmono3, errPr3)
	bK4 = coeffit_Kaiser(j, fcc, kstop,Pmm, biasF4, k, Pmono4, errPr4)
	#~ cn1 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt'
	#~ with open(cn1, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bK1[0][0],bK2[0][0],bK3[0][0],bK4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Scoccimaro'
	bsco1 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF1, k, Pmono1, errPr1)
	bsco2 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF2, k, Pmono2, errPr2)
	bsco3 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF3, k, Pmono3, errPr3)
	bsco4 = coeffit_Scocci(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF4, k, Pmono4, errPr4)
	#~ cn2 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt'
	#~ with open(cn2, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (bsco1[0][0],bsco2[0][0],bsco3[0][0],bsco4[0][0]))
	#~ fid_file.close()

	#----------------------------------------------------------------------------------------
	print 'Tns'
	btns1 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF1, k, Pmono1, errPr1, AB2_1, AB4_1, AB6_1, AB8_1)
	btns2 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF2, k, Pmono2, errPr2, AB2_2, AB4_2, AB6_2, AB8_2)
	btns3 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF3, k, Pmono3, errPr3, AB2_3, AB4_3, AB6_3, AB8_3)
	btns4 = coeffit_TNS(j, fcc, kstop,Pmm,Pmod_dt, Pmod_tt, biasF4, k, Pmono4, errPr4, AB2_4, AB4_4, AB6_4, AB8_4)
	#~ cn3 = '/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt'
	#~ with open(cn3, 'w') as fid_file:
		#~ fid_file.write('%.8g %.8g %.8g %.8g\n' % (btns1[0][0],btns2[0][0],btns3[0][0],btns4[0][0]))
	#~ fid_file.close()

	#~ #----------------------------------------------------------------------------------------
	print 'eTns'

	betns1 = coeffit_eTNS(j, fcc, kstop, b1pt3[0]*sca1, b2pt3[0]*sca1, bspt3[0]*sca1, b3pt3[0]*sca1, Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono1, errPr1, AB2bis_1, AB4bis_1,\
	AB6bis_1, AB8bis_1, Npt3[0])
	betns2 = coeffit_eTNS(j, fcc, kstop, b1pt3[1]*sca2, b2pt3[1]*sca2, bspt3[1]*sca2, b3pt3[1]*sca2, Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono2, errPr2, AB2bis_2, AB4bis_2,\
	AB6bis_2, AB8bis_2, Npt3[1])
	betns3 = coeffit_eTNS(j, fcc, kstop, b1pt3[2]*sca3, b2pt3[2]*sca3, bspt3[2]*sca3, b3pt3[2]*sca3, Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono3, errPr3, AB2bis_3, AB4bis_3,\
	AB6bis_3, AB8bis_3, Npt3[2])
	betns4 = coeffit_eTNS(j, fcc, kstop, b1pt3[3]*sca4, b2pt3[3]*sca4, bspt3[3]*sca4, b3pt3[3]*sca4, Pmm, Pmod_dt, Pmod_tt,\
	A, B, C, D, E, F, G, H, k, Pmono4, errPr4, AB2bis_4, AB4bis_4,\
	AB6bis_4, AB8bis_4, Npt3[3])
	
	
##########################################################################
##########################################################################

	#### compute the different power spectra given the mcmc results
	bK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')
	#~ BK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispkaibis_z='+str(z[j])+'.txt')
	#~ biasK = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispkai_z='+str(z[j])+'.txt')

	def kaips(b,sigma):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return Pmm*(b**2*coeffA +  2/3.*b*fcc*coeffB + 1/5.*fcc**2*coeffC)
		#~ return Pmm*(b**2 +  2/3.*b*fcc + 1/5.*fcc**2)
		
	#~ kai1 = kaips(biasF1, bK[0])
	#~ kai2 = kaips(biasF2, bK[1])
	#~ kai3 = kaips(biasF3, bK[2])
	#~ kai4 = kaips(biasF4, bK[3])
	kai1 = kaips(biasF1, bK1)
	kai2 = kaips(biasF2, bK2)
	kai3 = kaips(biasF3, bK3)
	kai4 = kaips(biasF4, bK4)
	#*****************************
	#~ k1 = kaips(bF1, BK[0])
	#~ k2 = kaips(bF2, BK[1])
	#~ k3 = kaips(bF3, BK[2])
	#~ k4 = kaips(bF4, BK[3])

	#~ k1ter = kaips(bF1, biasK[0])
	#~ k2ter = kaips(bF2, biasK[1])
	#~ k3ter = kaips(bF3, biasK[2])
	#~ k4ter = kaips(bF4, biasK[3])

	#---------------------------------------------------------------------------------------
	bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	#~ Bsco = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispscobis_z='+str(z[j])+'.txt')
	#~ bs = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispsco_z='+str(z[j])+'.txt')
	def scops(b,sigma):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		return b**2*Pmm*coeffA + 2/3.*b*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt

	#~ sco1 = scops(biasF1, bsco[0])
	#~ sco2 = scops(biasF2, bsco[1])
	#~ sco3 = scops(biasF3, bsco[2])
	#~ sco4 = scops(biasF4, bsco[3])
	sco1 = scops(biasF1, bsco1)
	sco2 = scops(biasF2, bsco2)
	sco3 = scops(biasF3, bsco3)
	sco4 = scops(biasF4, bsco4)
	
	#~ s1 = scops(bF1, Bsco[0])
	#~ s2 = scops(bF2, Bsco[1])
	#~ s3 = scops(bF3, Bsco[2])
	#~ s4 = scops(bF4, Bsco[3])
	
	#~ s1ter = scops(bF1, bs[0])
	#~ s2ter = scops(bF2, bs[1])
	#~ s3ter = scops(bF3, bs[2])
	#~ s4ter = scops(bF4, bs[3])
	#~ #---------------------------------------------------------------------------------------
	btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	#~ Btns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdisptnsbis_z='+str(z[j])+'.txt')
	#~ bt = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdisptns_z='+str(z[j])+'.txt')
	def tnsps(b,sigma, AB2, AB4, AB6, AB8):
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		return b**2*Pmm*coeffA + 2/3.*b*fcc*Pmod_dt*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
		+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return b**2*Pmm + 2/3.*b*fcc*Pmod_dt + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8)

	#~ tns1 = tnsps(biasF1,btns[0], AB2_1, AB4_1, AB6_1, AB8_1)
	#~ tns2 = tnsps(biasF2,btns[1], AB2_2, AB4_2, AB6_2, AB8_2)
	#~ tns3 = tnsps(biasF3,btns[2], AB2_3, AB4_3, AB6_3, AB8_3)
	#~ tns4 = tnsps(biasF4,btns[3], AB2_4, AB4_4, AB6_4, AB8_4)
	tns1 = tnsps(biasF1,btns1, AB2_1, AB4_1, AB6_1, AB8_1)
	tns2 = tnsps(biasF2,btns2, AB2_2, AB4_2, AB6_2, AB8_2)
	tns3 = tnsps(biasF3,btns3, AB2_3, AB4_3, AB6_3, AB8_3)
	tns4 = tnsps(biasF4,btns4, AB2_4, AB4_4, AB6_4, AB8_4)
	
	#~ t1 = tnsps(bF1,Btns[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2 = tnsps(bF2,Btns[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3 = tnsps(bF3,Btns[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4 = tnsps(bF4,Btns[3], ab2_4, ab4_4, ab6_4, ab8_4)
	
	#~ t1ter = tnsps(bF1,bt[0], ab2_1, ab4_1, ab6_1, ab8_1)
	#~ t2ter = tnsps(bF2,bt[1], ab2_2, ab4_2, ab6_2, ab8_2)
	#~ t3ter = tnsps(bF3,bt[2], ab2_3, ab4_3, ab6_3, ab8_3)
	#~ t4ter = tnsps(bF4,bt[3], ab2_4, ab4_4, ab6_4, ab8_4)
	#-------------------------------------------------------------------
	betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	#~ Betns = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/'+str(Mnu)+'eV/v_disp/case'+str(case)+'/vdispetnsbis_z='+str(z[j])+'.txt')
	#~ be = np.loadtxt('/home/david/codes/montepython_public/montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/v_disp/case'+str(case)+'/vdispetns_z='+str(z[j])+'.txt')
	def etnsps(b1,b2,bs,b3nl,sigma, AB2, AB4, AB6, AB8, N, sca=None):
		PsptD1z = b1**2*Pmm + b1*b2*A+ 1/4.*b2**2*B+ b1*bs*C+ 1/2.*b2*bs*D+ 1/4.*bs**2*E+ 2*b1*b3nl*F + N
		PsptT = b1* Pmod_dt+ b2*G+ bs*H + b3nl*F
		#~ kappa = k*sigma
		#~ coeffA = math.sqrt(math.pi)/2. * erf(kappa)/kappa
		#~ coeffB = 3./2./kappa**2*(coeffA - np.exp(-kappa**2))
		#~ coeffC = 5./2./kappa**2*(coeffB - np.exp(-kappa**2))
		#~ coeffD = 7./2./kappa**2*(coeffC - np.exp(-kappa**2))
		#~ coeffE = 9./2./kappa**2*(coeffD - np.exp(-kappa**2))
		kappa = k*sigma*fcc*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
		if sca:
			return  PsptD1z*coeffA*sca**2 + 2/3.*fcc*PsptT*coeffB*sca + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		else:
			return  PsptD1z*coeffA + 2/3.*fcc*PsptT*coeffB + 1/5.*fcc**2*Pmod_tt*coeffC \
			+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
		#~ return  PsptD1z + 2/3.*fcc*PsptT + 1/5.*fcc**2*Pmod_tt \
		#~ + (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
		
	
	#~ etns1 = etnsps(b1pt3[0], b2pt3[0], bspt3[0], b3pt3[0], betns[0], AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1)  
	#~ etns2 = etnsps(b1pt3[1], b2pt3[1], bspt3[1], b3pt3[1], betns[1], AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2)  
	#~ etns3 = etnsps(b1pt3[2], b2pt3[2], bspt3[2], b3pt3[2], betns[2], AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3)  
	#~ etns4 = etnsps(b1pt3[3], b2pt3[3], bspt3[3], b3pt3[3], betns[3], AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4) 
	etns1 = etnsps(b1pt3[0]*sca1, b2pt3[0]*sca1, bspt3[0]*sca1, b3pt3[0]*sca1, betns1, AB2bis_1, AB4bis_1, AB6bis_1, AB8bis_1, Npt3[0])  
	etns2 = etnsps(b1pt3[1]*sca2, b2pt3[1]*sca2, bspt3[1]*sca2, b3pt3[1]*sca2, betns2, AB2bis_2, AB4bis_2, AB6bis_2, AB8bis_2, Npt3[1])  
	etns3 = etnsps(b1pt3[2]*sca3, b2pt3[2]*sca3, bspt3[2]*sca3, b3pt3[2]*sca3, betns3, AB2bis_3, AB4bis_3, AB6bis_3, AB8bis_3, Npt3[2])  
	etns4 = etnsps(b1pt3[3]*sca4, b2pt3[3]*sca4, bspt3[3]*sca4, b3pt3[3]*sca4, betns4, AB2bis_4, AB4bis_4, AB6bis_4, AB8bis_4, Npt3[3]) 
	 
	#~ e1 = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], Betns[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1,(Bias_eff_t1/Bias_eff0_t1))  
	#~ e2 = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], Betns[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2,(Bias_eff_t2/Bias_eff0_t2))  
	#~ e3 = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], Betns[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3,(Bias_eff_t3/Bias_eff0_t3))  
	#~ e4 = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], Betns[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4,(Bias_eff_t4/Bias_eff0_t4))
	  
	#~ e1ter = etnsps(m1pt3[0], m1pt3[1], m1pt3[2], m1pt3[3], be[0], ab2bis_1, ab4bis_1, ab6bis_1, ab8bis_1)  
	#~ e2ter = etnsps(m2pt3[0], m2pt3[1], m2pt3[2], m2pt3[3], be[1], ab2bis_2, ab4bis_2, ab6bis_2, ab8bis_2)  
	#~ e3ter = etnsps(m3pt3[0], m3pt3[1], m3pt3[2], m3pt3[3], be[2], ab2bis_3, ab4bis_3, ab6bis_3, ab8bis_3)  
	#~ e4ter = etnsps(m4pt3[0], m4pt3[1], m4pt3[2], m4pt3[3], be[3], ab2bis_4, ab4bis_4, ab6bis_4, ab8bis_4)  

	#~ return kai1, kai2, kai3, kai4, 
	#~ return sco1, sco2, sco3, sco4
	#~ return tns1, tns2, tns3, tns4
	#~ return etns1, etns2, etns3, etns4
	return kai1, kai2, kai3, kai4, sco1, sco2, sco3, sco4, tns1, tns2, tns3, tns4, etns1, etns2, etns3, etns4
	
################################################################################################################
################################################################################################################
