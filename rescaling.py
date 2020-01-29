	def rescaling(self):
		
		dat_file_path = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/0.0eV/large_scale/'\
		'LS_z='+str(self.z)+'_.txt')
		with open(dat_file_path,'r') as f: 
				line = f.readline()   
				bcc_LS000 = float(line.split()[self.mbin])
				
		
		nu_masses = [0.03, 0.06, 0.1, 0.13, 0.15, 0.3, 0.45, 0.6]
		bcc_massive = np.zeros((len(nu_masses)))
		for count, mn in enumerate(nu_masses):
			dat_file_path2 = os.path.join(self.data_directory, 'montepython/likelihoods/BE_HaPPy/coefficients/other neutrinos masses/'+\
			str(mn)+'eV/LS_z='+str(self.z)+'_.txt')
			with open(dat_file_path2,'r') as f2: 
				line2 = f2.readline()   
				bcc_massive[count] = float(line2.split()[self.mbin])
				
		# interpolate on z
		bcc_final = np.interp(self.Mnu, nu_masses, bcc_massive)/bcc_LS000

		return bcc_final
