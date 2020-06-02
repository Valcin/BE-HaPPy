### code written by David Valcin
### for calibration details see arXiv:1901.06045

import numpy as np
import os
import math
import numpy as np
import sys
sys.path.append('/home/david/codes/FAST-PT')
import myFASTPT as FPT
dir_path = os.path.dirname(os.path.realpath(__file__))

def red_ps(mbin, bmodel, karray, z, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Plin, Pmod_dt,
	Pmod_tt, alpha, fog, rsd, red, kcase, sigma_v = []):

	####################################################################
	###### fit the Finger of God effect
	####################################################################
	#~ kbis = np.logspace(np.log10(np.min(kcamb)), np.log10(np.max(kcamb)), 200)
	#~ Plinbis = np.interp(kbis, k, Plin)

	#### compute tns coefficeints given mcmc results
	# set the parameters for the power spectrum window and
	# Fourier coefficient window 
	#P_window=np.array([.2,.2])  
	#~ C_window=0.95
	
	#~ for j in range(len(red)):
		#~ camb = np.loadtxt('/home/david/codes/Paco/data2/0.0eV/CAMB/Pk_cb_z='+str(red[j])+'00.txt')
		#~ kcamb = camb[:,0]
		#~ Pcamb = camb[:,1]
		#~ # padding length 
		#~ nu=-2; n_pad=len(kcamb)
		#~ n_pad=int(0.5*len(kcamb))
		#~ to_do=['all']
						
		#~ # initialize the FASTPT class 
		#~ # including extrapolation to higher and lower k  
		#~ # time the operation

		#~ fastpt=FPT.FASTPT(kcamb,to_do=to_do,n_pad=n_pad, verbose=True) 

		
		#~ AB2_1,AB4_1,AB6_1,AB8_1 = fastpt.RSD_ABsum_components(Pcamb,fz,b1 ,C_window=C_window)

		
		#~ dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV/TNS_coeff/'\
		#~ 'AB_'+str(mbin)+'_'+str(bmodel)+'_z='+str(red[j])+'.txt')	
		#~ with open(dat_file_path, 'w') as fid_file:
			#~ for index_k in xrange(len(kcamb)):
				#~ fid_file.write('%.8g %.8g %.8g %.8g %.8g\n' % (kcamb[index_k], AB2_1[index_k],AB4_1[index_k],AB6_1[index_k],AB8_1[index_k]))

	#~ kill


	AB2_temp = np.zeros((400,4))
	AB4_temp = np.zeros((400,4))
	AB6_temp = np.zeros((400,4))
	AB8_temp = np.zeros((400,4))
	AB2_temp2 = np.zeros(400)
	AB4_temp2 = np.zeros(400)
	AB6_temp2 = np.zeros(400)
	AB8_temp2 = np.zeros(400)
	
	for j in range(len(red)):
		dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV/TNS_coeff/'\
		'AB_'+str(mbin)+'_'+str(bmodel)+'_z='+str(red[j])+'.txt')	
		kab = np.loadtxt(dat_file_path)[:,0]
		lkab = len(kab)

		with open(dat_file_path,'r') as f: 
			line = f.readline()
			for index_k in range(lkab):
				AB2_temp[index_k,j] = float(line.split()[0])
				AB4_temp[index_k,j] = float(line.split()[1])
				AB6_temp[index_k,j] = float(line.split()[2])
				AB8_temp[index_k,j] = float(line.split()[3])
				line = f.readline()
	for k in range(lkab): # interp on redshift
		AB2_temp2[k] = np.interp(z,red,AB2_temp[k,:])
		AB4_temp2[k] = np.interp(z,red,AB4_temp[k,:])
		AB6_temp2[k] = np.interp(z,red,AB6_temp[k,:])
		AB8_temp2[k] = np.interp(z,red,AB8_temp[k,:])
		
	# interp on k array
	AB2 = np.interp(karray, kab, AB2_temp2) 
	AB4 = np.interp(karray, kab, AB4_temp2) 
	AB6 = np.interp(karray, kab, AB6_temp2)
	AB8 = np.interp(karray, kab, AB8_temp2) 
		
		
		
	if fog == 1:
		if rsd == 1:
			mname = 'kai'
		elif rsd == 2:
			mname = 'sco'
		elif rsd == 3 and bmodel == 1:
			mname = 'tns'
		elif rsd == 3 and bmodel == 2:
			mname = 'tns'
		elif rsd == 3 and bmodel == 3:
			mname = 'etns'

		if sigma_v == []:
			sigma_temp = np.zeros(4)
			for j in range(len(red)):
				#~ for count,iz in enumerate(red):
				dat_path = os.path.join(dir_path, 'coefficients/0.0eV'\
				'/v_disp/case'+str(kcase)+'/vdisp'+mname+'_z='+str(red[j])+'.txt')
				with open(dat_path,'r') as f:
					for i, line in enumerate(f):
							sigma_temp[j] = float(line.split()[mbin])
			sigma_v = np.interp(z,red,sigma_temp)
							

		# compute the multipole expansion coefficients
		kappa = np.zeros((len(karray)))
		coeffA = np.zeros((len(karray)))
		coeffB = np.zeros((len(karray)))
		coeffC = np.zeros((len(karray)))
		coeffD = np.zeros((len(karray)))
		coeffE = np.zeros((len(karray)))

		kappa = karray*sigma_v*fz*Dz
		coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
		coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
		coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
		coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
		coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))

		if rsd == 1:
			if bmodel == 1:
				b = b1 * alpha
				Pred = Plin*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = Plin*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
				
			elif bmodel == 3:
				raise ValueError('Sorry combination not available')
				
		elif rsd == 2:
			if bmodel == 1:
				b = b1 * alpha
				Pred = Plin*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = Plin*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
				
			elif bmodel == 3:
				raise ValueError('Sorry combination not available')

		elif rsd == 3:
			if bmodel == 1:
				b = b1 * alpha
				Pred = b**2*Plin*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
				+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = b**2*Plin*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
				+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
				
			elif bmodel == 3:
				# here b3 == bs and b4 == b3nl
				PhhDD = (b1**2*Plin + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
				2*b1*b4*F )* alpha**2
				PhhDT = (b1* Pmod_dt + b2*G + b3*H + b4*F) * alpha
				Pred = PhhDD*coeffA  + 2/3.*fz*PhhDT*coeffB + 1/5.*fz**2*Pmod_tt*coeffC + 1/3.*AB2*coeffB \
				+ 1/5.*AB4*coeffC + 1/7.*AB6*coeffD + 1/9.*AB8*coeffE 
					
	#------------
	else:
		if rsd == 1:
			if bmodel == 1:
				b = b1* alpha
				Pred = Plin*(b**2 + 2/3.*b*fz + 1/5.*fz**2) 
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = Plin*(b**2 + 2/3.*b*fz + 1/5.*fz**2 )
				
			elif bmodel == 3:
				raise ValueError('Sorry combination not available')
				
		elif rsd == 2:
			if bmodel == 1:
				b = b1* alpha
				Pred = Plin*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = Plin*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
				
			elif bmodel == 3:
				raise ValueError('Sorry combination not available')

		elif rsd == 3:
			if bmodel == 1:
				b = b1* alpha
				Pred = b**2*Plin + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
				+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
				
			elif bmodel == 2:
				b = (b1 + b2*(karray**2) + b3*(karray**3) + b4*(karray**4))*alpha
				Pred = b**2*Plin + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
				+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
				
			elif bmodel == 3:
				PhhDD = (b1**2*Plin + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
				2*b1*b4*F )*alpha
				PhhDT = (b1* Pmod_dt + b2*G + b3*H + b4*F )*alpha
				Pred = PhhDD  + 2/3.*fz*PhhDT + 1/5.*fz**2*Pmod_tt + 1/3.*AB2 \
				+ 1/5.*AB4 + 1/7.*AB6 + 1/9.*AB8 
		
	return Pred
