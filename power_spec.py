	def red_ps(self,data, kbound, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, Pmod_dt, Pmod_tt, alpha):
			
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/TNS_coeff/'\
		'AB_'+str(self.mbin)+'_'+str(self.bmodel)+'_z='+str(self.z)+'.txt')	
		AB2 = np.zeros((len(kbound)))
		AB4 = np.zeros((len(kbound)))
		AB6 = np.zeros((len(kbound)))
		AB8 = np.zeros((len(kbound)))
		with open(dat_file_path,'r') as f: 
			line = f.readline()
			for index_k in range(len(kbound)):
				AB2[index_k] = float(line.split()[0])
				AB4[index_k] = float(line.split()[1])
				AB6[index_k] = float(line.split()[2])
				AB8[index_k] = float(line.split()[3])
				line = f.readline()
			
		if self.fog == 1:
			# compute the multipole expansion coefficients
			kappa = np.zeros((len(kbound)))
			coeffA = np.zeros((len(kbound)))
			coeffB = np.zeros((len(kbound)))
			coeffC = np.zeros((len(kbound)))
			coeffD = np.zeros((len(kbound)))
			coeffE = np.zeros((len(kbound)))
			
			sigma_v = (data.mcmc_parameters['sigma_v']['current'] *
			data.mcmc_parameters['sigma_v']['scale'])

			kappa = kbound*sigma_v*fz*Dz
			coeffA = np.arctan(kappa/math.sqrt(2))/(math.sqrt(2)*kappa) + 1/(2+kappa**2)
			coeffB = 6/kappa**2*(coeffA - 2/(2+kappa**2))
			coeffC = -10/kappa**2*(coeffB - 2/(2+kappa**2))
			coeffD = -2/3./kappa**2*(coeffC - 2/(2+kappa**2))
			coeffE = -4/10./kappa**2*(7.*coeffD - 2/(2+kappa**2))
	
			if self.rsd == 1:
				if self.bmodel == 1:
					b = b1 * alpha
					Pred = Pmod_dd*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
					
				elif self.bmodel == 2:
					b = (b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4))*alpha
					Pred = Pmod_dd*(b**2*coeffA + 2/3.*b*fz*coeffB + 1/5.*fz**2*coeffC) 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
					
			elif self.rsd == 2:
				if self.bmodel == 1:
					b = b1 * alpha
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
					
				elif self.bmodel == 2:
					b = (b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4))*alpha
					Pred = Pmod_dd*b**2*coeffA + 2/3.*b*fz*coeffB*Pmod_dt + 1/5.*fz**2*coeffC*Pmod_tt 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
	
			elif self.rsd == 3:
				if self.bmodel == 1:
					b = b1
					Pred = b**2*Pmod_dd*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
					+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = b**2*Pmod_dd*coeffA + 2/3.*b*fz*Pmod_dt*coeffB + 1/5.*fz**2*Pmod_tt*coeffC \
					+ (1/3.*AB2*coeffB+ 1/5.*AB4*coeffC+ 1/7.*AB6*coeffD+ 1/9.*AB8*coeffE)
					
				elif self.bmodel == 3:
					# here b3 == bs and b4 == b3nl
					PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
					2*b1*b4*F 
					PhhDT = b1* Pmod_dt + b2*G + b3*H + b4*F 
					Pred = PhhDD*coeffA  + 2/3.*fz*PhhDT*coeffB + 1/5.*fz**2*Pmod_tt*coeffC + 1/3.*AB2*coeffB \
					+ 1/5.*AB4*coeffC + 1/7.*AB6*coeffD + 1/9.*AB8*coeffE 
						
		#------------
		else:
			if self.rsd == 1:
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*(b**2 + 2/3.*b*fz + 1/5.*fz**2) 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = Pmod_dd*(b**2 + 2/3.*b*fz + 1/5.*fz**2 )
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
					
			elif self.rsd == 2:
				if self.bmodel == 1:
					b = b1
					Pred = Pmod_dd*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = Pmod_dd*b**2 + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt 
					
				elif self.bmodel == 3:
					raise ValueError('Sorry combination not available')
	
			elif self.rsd == 3:
				if self.bmodel == 1:
					b = b1
					Pred = b**2*Pmod_dd + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
					+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
					
				elif self.bmodel == 2:
					b = b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4) 
					Pred = b**2*Pmod_dd + 2/3.*b*fz*Pmod_dt + 1/5.*fz**2*Pmod_tt \
					+ (1/3.*AB2+ 1/5.*AB4+ 1/7.*AB6+ 1/9.*AB8) 
					
				elif self.bmodel == 3:
					# here b3 == bs and b4 == b3nl
					PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
					2*b1*b4*F 
					PhhDT = b1* Pmod_dt + b2*G + b3*H + b4*F 
					Pred = PhhDD  + 2/3.*fz*PhhDT + 1/5.*fz**2*Pmod_tt + 1/3.*AB2 \
					+ 1/5.*AB4 + 1/7.*AB6 + 1/9.*AB8 
			
		return Pred
