
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
		cpart1 = '/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='
		cpart2 = '/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='
		cpart3 = '/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='
		cpart4 = '/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh1_realisation_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh2_realisation_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh3_realisation_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/Phh4_realisation_z='+str(2.0)+'.txt')
	elif mv == 0.15:
		cpart1 = '/home/david/codes/Paco/data2/0.15eV/Phh1_realisation_0.15_z='
		cpart2 = '/home/david/codes/Paco/data2/0.15eV/Phh2_realisation_0.15_z='
		cpart3 = '/home/david/codes/Paco/data2/0.15eV/Phh3_realisation_0.15_z='
		cpart4 = '/home/david/codes/Paco/data2/0.15eV/Phh4_realisation_0.15_z='
		#~ d1 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh1_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d2 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh2_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d3 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh3_realisation_0.15_z='+str(2.0)+'.txt')
		#~ d4 = np.loadtxt('/home/david/codes/Paco/data2/0.15eV/Phh4_realisation_0.15_z='+str(2.0)+'.txt')
	
	####################################################################
	d1 = np.loadtxt(cpart1+str(0.0)+'.txt')
	k = d1[:,19]
	P1 = np.zeros((len(k),l2))
	P2 = np.zeros((len(k),l2))
	P3 = np.zeros((len(k),l2))
	P4 = np.zeros((len(k),l2))
	
	for j, jz in enumerate(red2):
		d1 = np.loadtxt(cpart1+str(jz)+'.txt')
		d2 = np.loadtxt(cpart2+str(jz)+'.txt')
		d3 = np.loadtxt(cpart3+str(jz)+'.txt')
		d4 = np.loadtxt(cpart4+str(jz)+'.txt')
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
	
		P1[:,j] = PH1
		P2[:,j] = PH2
		P3[:,j] = PH3
		P4[:,j] = PH4
	
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
	#~ plt.plot(kclass,P1bis [:,3], color='k')
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
