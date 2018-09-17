#### Link to Alvise article 1704
#### Link to VILLAESCUSA -NAVARRO article 1708
#### Link to Hernandez article 1608
#### Link to Linder article 1211

from classy import Class
from matplotlib.colors import LogNorm
from cimport import classy_import
from bcoeff import bcoeff
from ls_coeff import lscoeff
from pt_coeff import ptcoeff
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
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	####################################################################
	#### select the k mode according to the ones in Raccanelli et al. 2017
	k_article = [0.35,0.2,0.15]
	
	if case == 1:
		kmax = k_article[0]
	elif case == 2:
		kmax = k_article[1]
	elif case == 3:
		kmax = k_article[2]
		
	####################################################################
	#### import classy results
	redshift, mv, h, d_tot, T_cb, pk, f, D = classy_import(self,cosmo, data, red2 )

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
	if model == 'lin':
		# get the fitted ps from p_halo
		k, P_halo = Halo(self, cosmo, data, 'lin', case, Massbins)
		
		# tinker effective bias
		bcc = lscoeff(self,data, mv,Massbins)[1]
					
		# compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.
		bmm = np.zeros((len(k),l2, len(Massbins)), 'float64')
		for iz in xrange(l2):
			for count,j in enumerate(Massbins):
				bmm[:,iz, count] = bcc[iz,count] * (T_cb[:,iz]/d_tot[:,iz])
				
		print np.shape(bmm)

		# Compute the halo Power spectrum in real space
		if RSD == 0:
			print 'you chose Kaiser model'
			dim = np.shape(P_halo)
			Pred = np.zeros(dim)
			for iz in xrange(znumber):
				for count,j in enumerate(Massbins):
					ind2 = mbins.index(j)
					Pred[:,iz,count] = P_halo[:,iz,count] + 2/3.*bmm[:,iz, count]*f[iz] + 1/5.*f[iz]**2

		return k, Pred
	#---------------------------------------------------------------
	#~ if fog:
		#~ # compute the halo ps in redshift space
		#~ Phh = np.zeros((len(kbis),znumber,len(Massbins)))
		#~ for iz in xrange(znumber):
			#~ for count,j in enumerate(Massbins):
				#~ ind2 = mbins.index(j)
				#~ # interpolate on the new scale array
				#~ pk_lin2[:,iz] = np.interp(kbis,k, pk_lin[:,iz])
				#~ Pdd2[:,iz,count] = np.interp(kbis,k, Pdd[:,iz,count])
				#~ Pdt2[:,iz,count] = np.interp(kbis,k, Pdt[:,iz,count])
				#~ Ptt2[:,iz] = np.interp(kbis,k, Ptt[:,iz])
				#~ AB2,AB4,AB6,AB8 = fastpt2.RSD_ABsum_components(pk_lin2[:,iz],f[iz], b1[iz,count],C_window=C_window) #tns coeff
				#~ Phh[:,iz,count] = Pdd2[:,iz,count]*coeffA + 2/3.*f[iz]*Pdt2[:,iz,count]*coeffB +\
				#~ 1/5.*f[iz]**2*Ptt2[:,iz]*coeffC + 1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE
	#~ else:
		#~ Phh = np.zeros((len(kbis),znumber,len(Massbins)))
		#~ for iz in xrange(znumber):
			#~ for count,j in enumerate(Massbins):
				#~ ind2 = mbins.index(j)
				#~ # interpolate on the new scale array
				#~ pk_lin2[:,iz] = np.interp(kbis,k, pk_lin[:,iz])
				#~ Pdd2[:,iz,count] = np.interp(kbis,k, Pdd[:,iz,count])
				#~ Pdt2[:,iz,count] = np.interp(kbis,k, Pdt[:,iz,count])
				#~ Ptt2[:,iz] = np.interp(kbis,k, Ptt[:,iz])
				#~ AB2,AB4,AB6,AB8 = fastpt.RSD_ABsum_components(pk_lin2[:,iz],f[iz], b1[iz,count],C_window=C_window) #tns coeff
				#~ Phh[:,iz,count] = Pdd2[:,iz,count] + 2/3.*f[iz]*Pdt2[:,iz,count] +\
				#~ 1/5.*f[iz]**2*Ptt2[:,iz] + 1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6 + 1/9.*AB8
					
		
		#---------------------------------------------------
		#--------- halo redshift space -------------------------
		#---------------------------------------------------
		#~ d1a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_0_z='+str(2.0)+'.txt')
		#~ d1b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_0_z='+str(2.0)+'.txt')
		#~ d1c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_0_z='+str(2.0)+'.txt')
		#~ d1d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_0_z='+str(2.0)+'.txt')
		#~ d2a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_1_z='+str(2.0)+'.txt')
		#~ d2b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_1_z='+str(2.0)+'.txt')
		#~ d2c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_1_z='+str(2.0)+'.txt')
		#~ d2d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_1_z='+str(2.0)+'.txt')
		#~ d3a = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_red_axis_2_z='+str(2.0)+'.txt')
		#~ d3b = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_red_axis_2_z='+str(2.0)+'.txt')
		#~ d3c = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_red_axis_2_z='+str(2.0)+'.txt')
		#~ d3d = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_red_axis_2_z='+str(2.0)+'.txt')


		#~ kx1a = d1a[:,19]
		#~ Px1a = np.zeros((len(kx1a),10))
		#~ Px1b = np.zeros((len(kx1a),10))
		#~ Px1c = np.zeros((len(kx1a),10))
		#~ Px1d = np.zeros((len(kx1a),10))
		#~ Pxshot1a = np.zeros((10))
		#~ Pxshot1b = np.zeros((10))
		#~ Pxshot1c = np.zeros((10))
		#~ Pxshot1d = np.zeros((10))
		#~ pnum1 = [0,2,4,6,8,10,12,14,16,18]
		#~ pnum2 = [1,3,5,7,9,11,13,15,17,20]
		#~ for i in xrange(0,10):
			#~ Px1a[:,i]= d1a[:,pnum1[i]]
			#~ Px1b[:,i]= d1b[:,pnum1[i]]
			#~ Px1c[:,i]= d1c[:,pnum1[i]]
			#~ Px1d[:,i]= d1d[:,pnum1[i]]
			#~ Pxshot1a[i]= d1a[0,pnum2[i]]
			#~ Pxshot1b[i]= d1b[0,pnum2[i]]
			#~ Pxshot1c[i]= d1c[0,pnum2[i]]
			#~ Pxshot1d[i]= d1d[0,pnum2[i]]
		#~ kx2a = d2a[:,19]
		#~ Px2a = np.zeros((len(kx2a),10))
		#~ Px2b = np.zeros((len(kx2a),10))
		#~ Px2c = np.zeros((len(kx2a),10))
		#~ Px2d = np.zeros((len(kx2a),10))
		#~ Pxshot2a = np.zeros((10))
		#~ Pxshot2b = np.zeros((10))
		#~ Pxshot2c = np.zeros((10))
		#~ Pxshot2d = np.zeros((10))
		#~ pnum1 = [0,2,4,6,8,10,12,14,16,18]
		#~ pnum2 = [1,3,5,7,9,11,13,15,17,20]
		#~ for i in xrange(0,10):
			#~ Px2a[:,i]= d2a[:,pnum1[i]]
			#~ Px2b[:,i]= d2b[:,pnum1[i]]
			#~ Px2c[:,i]= d2c[:,pnum1[i]]
			#~ Px2d[:,i]= d2d[:,pnum1[i]]
			#~ Pxshot2a[i]= d2a[0,pnum2[i]]
			#~ Pxshot2b[i]= d2b[0,pnum2[i]]
			#~ Pxshot2c[i]= d2c[0,pnum2[i]]
			#~ Pxshot2d[i]= d2d[0,pnum2[i]]
		#~ kx3a = d3a[:,19]
		#~ Px3a = np.zeros((len(kx3a),10))
		#~ Px3b = np.zeros((len(kx3a),10))
		#~ Px3c = np.zeros((len(kx3a),10))
		#~ Px3d = np.zeros((len(kx3a),10))
		#~ Pxshot3a = np.zeros((10))
		#~ Pxshot3b = np.zeros((10))
		#~ Pxshot3c = np.zeros((10))
		#~ Pxshot3d = np.zeros((10))
		#~ pnum1 = [0,2,4,6,8,10,12,14,16,18]
		#~ pnum2 = [1,3,5,7,9,11,13,15,17,20]
		#~ for i in xrange(0,10):
			#~ Px3a[:,i]= d3a[:,pnum1[i]]
			#~ Px3b[:,i]= d3b[:,pnum1[i]]
			#~ Px3c[:,i]= d3c[:,pnum1[i]]
			#~ Px3d[:,i]= d3d[:,pnum1[i]]
			#~ Pxshot3a[i]= d3a[0,pnum2[i]]
			#~ Pxshot3b[i]= d3b[0,pnum2[i]]
			#~ Pxshot3c[i]= d3c[0,pnum2[i]]
			#~ Pxshot3d[i]= d3d[0,pnum2[i]]
			
		#~ for i in xrange(0,10):
			#~ Px1a[:,i] = Px1a[:,i]-Pxshot1a[i]
			#~ Px1b[:,i] = Px1b[:,i]-Pxshot1b[i]
			#~ Px1c[:,i] = Px1c[:,i]-Pxshot1c[i]
			#~ Px1d[:,i] = Px1d[:,i]-Pxshot1d[i]
			#~ Px2a[:,i] = Px2a[:,i]-Pxshot2a[i]
			#~ Px2b[:,i] = Px2b[:,i]-Pxshot2b[i]
			#~ Px2c[:,i] = Px2c[:,i]-Pxshot2c[i]
			#~ Px2d[:,i] = Px2d[:,i]-Pxshot2d[i]
			#~ Px3a[:,i] = Px3a[:,i]-Pxshot3a[i]
			#~ Px3b[:,i] = Px3b[:,i]-Pxshot3b[i]
			#~ Px3c[:,i] = Px3c[:,i]-Pxshot3c[i]
			#~ Px3d[:,i] = Px3d[:,i]-Pxshot3d[i]
			
			#~ nul1a = np.where(Px1a[:,i] < 0)[0]
			#~ Px1a[nul1a,i] = 0
			#~ nul1b = np.where(Px1b[:,i] < 0)[0]
			#~ Px1b[nul1b,i] = 0
			#~ nul1c = np.where(Px1c[:,i] < 0)[0]
			#~ Px1c[nul1c,i] = 0
			#~ nul1d = np.where(Px1d[:,i] < 0)[0]
			#~ Px1d[nul1d,i] = 0
			#~ nul2a = np.where(Px2a[:,i] < 0)[0]
			#~ Px2a[nul2a,i] = 0
			#~ nul2b = np.where(Px2b[:,i] < 0)[0]
			#~ Px2b[nul2b,i] = 0
			#~ nul2c = np.where(Px2c[:,i] < 0)[0]
			#~ Px2c[nul2c,i] = 0
			#~ nul2d = np.where(Px2d[:,i] < 0)[0]
			#~ Px2d[nul2d,i] = 0
			#~ nul3a = np.where(Px3a[:,i] < 0)[0]
			#~ Px3a[nul3a,i] = 0
			#~ nul3b = np.where(Px3b[:,i] < 0)[0]
			#~ Px3b[nul3b,i] = 0
			#~ nul3c = np.where(Px3c[:,i] < 0)[0]
			#~ Px3c[nul3c,i] = 0
			#~ nul3d = np.where(Px3d[:,i] < 0)[0]
			#~ Px3d[nul3d,i] = 0

		#~ Pmono1temp = (Px1a + Px2a + Px3a)/3
		#~ Pmono2temp = (Px1b + Px2b + Px3b)/3
		#~ Pmono3temp = (Px1c + Px2c + Px3c)/3
		#~ Pmono4temp = (Px1d + Px2d + Px3d)/3


		#~ ### do the mean and std over quantitites ###
		
		#~ Pmono1 = np.mean(Pmono1temp[:,0:11], axis=1)
		#~ Pmono2 = np.mean(Pmono2temp[:,0:11], axis=1)
		#~ Pmono3 = np.mean(Pmono3temp[:,0:11], axis=1)
		#~ Pmono4 = np.mean(Pmono4temp[:,0:11], axis=1)

		#~ return kbis,Phh, kx1a, Pmono1, Pmono2, Pmono3, Pmono4

	



