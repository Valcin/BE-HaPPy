import numpy as np
import os
import math

dir_path = os.path.dirname(os.path.realpath(__file__))

def red_ps(mbin, bmodel, karray, z, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Plin, Pmod_dt,
	Pmod_tt, alpha, fog, rsd, red, kcase, sigma_v = []):


	dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV/TNS_coeff/'\
	'AB_'+str(mbin)+'_'+str(bmodel)+'_z='+str(z)+'.txt')	
	AB2 = np.zeros((len(karray)))
	AB4 = np.zeros((len(karray)))
	AB6 = np.zeros((len(karray)))
	AB8 = np.zeros((len(karray)))
	with open(dat_file_path,'r') as f: 
		line = f.readline()
		for index_k in range(len(karray)):
			AB2[index_k] = float(line.split()[0])
			AB4[index_k] = float(line.split()[1])
			AB6[index_k] = float(line.split()[2])
			AB8[index_k] = float(line.split()[3])
			line = f.readline()
		
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
			#~ for count,iz in enumerate(red):
			dat_path = os.path.join(dir_path, 'coefficients/0.0eV'\
			'/v_disp/case'+str(kcase)+'/vdisp'+mname+'_z='+str(z)+'.txt')
			with open(dat_path,'r') as f:
				for i, line in enumerate(f):
						sigma_v = float(line.split()[mbin])
							

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
