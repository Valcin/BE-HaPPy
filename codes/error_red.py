
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


def error(self,data, kclass, P_halo, redshift, mv, Massbins):
	
	# create the directory if not already done
	self.data_directory = data.path['root']
	directory = os.path.join(self.data_directory, 'BE_HaPPy/Errors')

	if not os.path.exists(directory):
		os.makedirs(directory)
	
	# remove older files before computing new ones
	fileList = os.listdir(directory)
	for fileName in fileList:
		os.remove(directory+"/"+fileName)
	
	
	kmin = np.min(kclass)
	kmax = np.max(kclass)
	print kmin, kmax
	####################################################################
	#### Store the selected redshifts in a array and deduce its length for the loops
	#### array manipulation because len() and size only work for znumber >1
	znumber = len(redshift)

	####################################################################
    #### Store the redshifts where bcc fit  and bcc Ls are available in arrays
	red2 = [0.0,0.5,1.0,2.0]
	l2= len(red2)
	
	####################################################################
	#### Store the mass bins available in an array
	mbins = ['M1', 'M2', 'M3', 'M4']
	
	if mv == 0.0:
		cpart1 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.0eV/Phh1_realisation_'
		cpart2 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.0eV/Phh2_realisation_'
		cpart3 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.0eV/Phh3_realisation_'
		cpart4 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.0eV/Phh4_realisation_'

	elif mv == 0.15:
		cpart1 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.15eV/Phh1_realisation_0.15_'
		cpart2 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.15eV/Phh2_realisation_0.15_'
		cpart3 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.15eV/Phh3_realisation_0.15_'
		cpart4 = '/home/david/codes/montepython_public/BE_HaPPy/coefficients/simu_ps/0.15eV/Phh4_realisation_0.15_'
	
	####################################################################
	d1 = np.loadtxt(cpart1+'red_axis_0_z='+str(0.0)+'.txt')
	k = d1[:,19]
	P1 = np.zeros((len(k),l2))
	P2 = np.zeros((len(k),l2))
	P3 = np.zeros((len(k),l2))
	P4 = np.zeros((len(k),l2))
	
	for j, jz in enumerate(red2):
		d1a = np.loadtxt(cpart1+'red_axis_0_z='+str(jz)+'.txt')
		d2a = np.loadtxt(cpart1+'red_axis_1_z='+str(jz)+'.txt')
		d3a = np.loadtxt(cpart1+'red_axis_2_z='+str(jz)+'.txt')
		#-----------------------------------------------------------------------------------------------------
		d1b = np.loadtxt(cpart2+'red_axis_0_z='+str(jz)+'.txt')
		d2b = np.loadtxt(cpart2+'red_axis_1_z='+str(jz)+'.txt')
		d3b = np.loadtxt(cpart2+'red_axis_2_z='+str(jz)+'.txt')
		#-----------------------------------------------------------------------------------------------------
		d1c = np.loadtxt(cpart3+'red_axis_0_z='+str(jz)+'.txt')
		d2c = np.loadtxt(cpart3+'red_axis_1_z='+str(jz)+'.txt')
		d3c = np.loadtxt(cpart3+'red_axis_2_z='+str(jz)+'.txt')
		#-----------------------------------------------------------------------------------------------------
		d1d = np.loadtxt(cpart4+'red_axis_0_z='+str(jz)+'.txt')
		d2d = np.loadtxt(cpart4+'red_axis_1_z='+str(jz)+'.txt')
		d3d = np.loadtxt(cpart4+'red_axis_2_z='+str(jz)+'.txt')
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
		

		P1[:,j] = Pmono1
		P2[:,j] = Pmono2
		P3[:,j] = Pmono3
		P4[:,j] = Pmono4


		
	####################################################################
	# interpolate on selected redshift
	Phh1bis = np.zeros((len(k),znumber))
	Phh2bis = np.zeros((len(k),znumber))
	Phh3bis = np.zeros((len(k),znumber))
	Phh4bis = np.zeros((len(k),znumber))
	for ik in xrange(len(k)):
		f1 = interp1d(red2, P1[ik,:], kind='cubic', fill_value='extrapolate')
		f2 = interp1d(red2, P2[ik,:], kind='cubic', fill_value='extrapolate')
		f3 = interp1d(red2, P3[ik,:], kind='cubic', fill_value='extrapolate')
		f4 = interp1d(red2, P4[ik,:], kind='cubic', fill_value='extrapolate')
		Phh1bis[ik,:] = f1(redshift)
		Phh2bis[ik,:] = f2(redshift)
		Phh3bis[ik,:] = f3(redshift)
		Phh4bis[ik,:] = f4(redshift)

	# interpolate on selected k
	P1bis = np.zeros((len(kclass),znumber))
	P2bis = np.zeros((len(kclass),znumber))
	P3bis = np.zeros((len(kclass),znumber))
	P4bis = np.zeros((len(kclass),znumber))
	for iz in xrange(znumber):
		f1 = interp1d(k, Phh1bis[:, iz], kind='cubic', fill_value='extrapolate')
		f2 = interp1d(k, Phh2bis[:, iz], kind='cubic', fill_value='extrapolate')
		f3 = interp1d(k, Phh3bis[:, iz], kind='cubic', fill_value='extrapolate')
		f4 = interp1d(k, Phh4bis[:, iz], kind='cubic', fill_value='extrapolate')
		P1bis[:, iz] = f1(kclass)
		P2bis[:, iz] = f2(kclass)
		P3bis[:, iz] = f3(kclass)
		P4bis[:, iz] = f4(kclass)
    

    
    
	### save the errors in txt format
	for count,j in enumerate(Massbins):
		ind2 = mbins.index(j)
		if ind2 == 0:
			err1 = np.zeros((len(kclass),znumber))
			for iz,rd in enumerate(redshift):
				for i in xrange(len(kclass)):
					err1[i,iz] = np.absolute(P_halo[i,iz,count] - P1bis[i,iz])/P1bis[i,iz]
				cname1 = os.path.join(self.data_directory, 'BE_HaPPy/Errors/Err1_'+str(mv)+'_z='+str(rd)+'_.txt')
				with open(cname1, 'w+') as fid_file:
					for index_k in xrange(len(kclass)):
						fid_file.write('%.8g %.8g %.8g\n' % ( kclass[index_k], err1[index_k, iz], err1[index_k, iz]*100. ))
					fid_file.close()
		elif ind2 == 1:
			err2 = np.zeros((len(kclass),znumber))
			for iz,rd in enumerate(redshift):
				for i in xrange(len(kclass)):
					err2[i,iz] = np.absolute(P_halo[i,iz,count] - P2bis[i,iz])/P2bis[i,iz]
				cname2 = os.path.join(self.data_directory, 'BE_HaPPy/Errors/Err2_'+str(mv)+'_z='+str(rd)+'_.txt')
				with open(cname2, 'w+') as fid_file:
					for index_k in xrange(len(kclass)):
						fid_file.write('%.8g %.8g %.8g\n' % ( kclass[index_k], err2[index_k, iz], err2[index_k, iz]*100. ))
					fid_file.close()
		elif ind2 == 2:
			err3 = np.zeros((len(kclass),znumber))
			for iz,rd in enumerate(redshift):
				for i in xrange(len(kclass)):
					err3[i,iz] = np.absolute(P_halo[i,iz,count] - P3bis[i,iz])/P3bis[i,iz]
				cname3 = os.path.join(self.data_directory, 'BE_HaPPy/Errors/Err3_'+str(mv)+'_z='+str(rd)+'_.txt')
				with open(cname3, 'w+') as fid_file:
					for index_k in xrange(len(kclass)):
						fid_file.write('%.8g %.8g %.8g\n' % ( kclass[index_k], err3[index_k, iz], err3[index_k, iz]*100. ))
					fid_file.close()
		elif ind2 == 3:
			err4 = np.zeros((len(kclass),znumber))
			for iz,rd in enumerate(redshift):
				for i in xrange(len(kclass)):
					err4[i,iz] = np.absolute(P_halo[i,iz,count] - P4bis[i,iz])/P4bis[i,iz]
				cname4 = os.path.join(self.data_directory, 'BE_HaPPy/Errors/Err4_'+str(mv)+'_z='+str(rd)+'_.txt')
				with open(cname4, 'w+') as fid_file:
					for index_k in xrange(len(kclass)):
						fid_file.write('%.8g %.8g %.8g\n' % ( kclass[index_k], err4[index_k, iz], err4[index_k, iz]*100. ))
					fid_file.close()
    
    
	####################################################################
	#### halo ps vs k, same redshift but different massbins
	
	#~ plt.figure()
	#~ plt.plot(kclass,P_halo[:,3,0])
	#~ plt.plot(kclass,P_halo[:,3,1])
	#~ plt.plot(kclass,P_halo[:,3,2])
	#~ plt.plot(kclass,P_halo[:,3,3])
	#~ plt.plot(kclass,P1bis[:,3], color='k')
	#~ plt.plot(kclass,P2bis[:,3], color='k')
	#~ plt.plot(kclass,P3bis[:,3], color='k')
	#~ plt.plot(kclass,P4bis[:,3], color='k')
	#~ plt.xscale('log')
	#~ plt.yscale('log')
	#~ plt.ylim(1e2,1e7)
	#~ plt.xlim(kmin,kmax)
	#~ plt.show()

	#### halo ps vs k, same massbin but different redshifts
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
