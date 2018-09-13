
from classy import Class
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
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


def error(kclass, P_halo,mv):
	
	if mv == 0.0:
		d1 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='+str(2.0)+'.txt')
		d2 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='+str(2.0)+'.txt')
		d3 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='+str(2.0)+'.txt')
		d4 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='+str(2.0)+'.txt')
	elif mv == 0.15:
		d1 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh1_realisation_0.15_z='+str(2.0)+'.txt')
		d2 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh2_realisation_0.15_z='+str(2.0)+'.txt')
		d3 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh3_realisation_0.15_z='+str(2.0)+'.txt')
		d4 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh4_realisation_0.15_z='+str(2.0)+'.txt')
	
	k = d1[:,19]
	Phh1 = np.zeros((len(k),10))
	Phh2 = np.zeros((len(k),10))
	Phh3 = np.zeros((len(k),10))
	Phh4 = np.zeros((len(k),10))
	Pshot1 = np.zeros((10))
	Pshot2 = np.zeros((10))
	Pshot3 = np.zeros((10))
	Pshot4 = np.zeros((10))
	pnum1 = [0,2,4,6,8,10,12,14,16,18]
	pnum2 = [1,3,5,7,9,11,13,15,17,20]
	for i in xrange(0,10):
		Phh1[:,i]= d1[:,pnum1[i]]
		Phh2[:,i]= d2[:,pnum1[i]]
		Phh3[:,i]= d3[:,pnum1[i]]
		Phh4[:,i]= d4[:,pnum1[i]]
		Pshot1[i]= d1[0,pnum2[i]]
		Pshot2[i]= d2[0,pnum2[i]]
		Pshot3[i]= d3[0,pnum2[i]]
		Pshot4[i]= d4[0,pnum2[i]]
		Phh1[:,i] = Phh1[:,i]-Pshot1[i]
		Phh2[:,i] = Phh2[:,i]-Pshot2[i]
		Phh3[:,i] = Phh3[:,i]-Pshot3[i]
		Phh4[:,i] = Phh4[:,i]-Pshot4[i]
		
	PH1 = np.mean(Phh1[:,0:11], axis=1)
	PH2 = np.mean(Phh2[:,0:11], axis=1)
	PH3 = np.mean(Phh3[:,0:11], axis=1)
	PH4 = np.mean(Phh4[:,0:11], axis=1)
	
	####################################################################
	#### halo ps vs k, same redshift but different massbins
	plt.figure()
	plt.plot(kclass,P_halo[:,3,0])
	plt.plot(kclass,P_halo[:,3,1])
	plt.plot(kclass,P_halo[:,3,2])
	plt.plot(kclass,P_halo[:,3,3])
	plt.plot(k,psimu1, color='k')
	plt.plot(k,psimu2, color='k')
	plt.plot(k,psimu3, color='k')
	plt.plot(k,psimu4, color='k')
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(1e2,1e7)
	plt.xlim(1e-3,0.2)
	plt.show()
	
	#### halo ps vs k, same massbin but different redshift
	#~ plt.figure()
	#~ plt.plot(kclass,P_halo[:,0,3])
	#~ plt.plot(kclass,P_halo[:,1,3])
	#~ plt.plot(kclass,P_halo[:,2,3])
	#~ plt.plot(kclass,P_halo[:,3,3])
	#~ plt.xscale('log')
	#~ plt.yscale('log')
	#~ plt.ylim(1e3,2e4)
	#~ plt.xlim(1e-3,0.2)
	#~ plt.show()
	
	return kclass, Phh, k, PH1, PH2, PH3, PH4
	
