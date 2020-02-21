import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

def bcoeff(mbin, bmodel, z, red, lred, kcase):
	b1 = np.zeros((lred))
	b2 = np.zeros((lred))
	b3 = np.zeros((lred))
	b4 = np.zeros((lred))
	b1_final = np.zeros(len(np.atleast_1d(z)))
	b2_final = np.zeros(len(np.atleast_1d(z)))
	b3_final = np.zeros(len(np.atleast_1d(z)))
	b4_final = np.zeros(len(np.atleast_1d(z)))
	
	
	
	if bmodel == 1:
		for count,iz in enumerate(red):
			dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV/large_scale/'\
			'LS_z='+str(iz)+'_.txt')
			with open(dat_file_path,'r') as f: 
				line = f.readline()   
				b1[count] = float(line.split()[mbin])
				
		# interpolate on z
		b1_final = np.interp(z, red, b1)
		b2_final = np.interp(z, red, b2)
		b3_final = np.interp(z, red, b3)
		b4_final = np.interp(z, red, b4)
		
		print(b1_final, b2_final, b3_final, b4_final)
		
		return b1_final, b2_final, b3_final, b4_final # here b2_final, b3_final, b4_final == 0

	elif bmodel == 2:
		for count,iz in enumerate(red):
			dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV'\
			'/case'+str(kcase)+'/coeff_pl_0.0_z='+str(iz)+'.txt')
			with open(dat_file_path,'r') as f:
				for i, line in enumerate(f):
					if i == mbin: 
						b1[count] = float(line.split()[0])
						b2[count] = float(line.split()[1])
						b3[count] = float(line.split()[2])
						b4[count] = float(line.split()[3])
						
		# interpolate on z
		b1_final = np.interp(z, red, b1)
		b2_final = np.interp(z, red, b2)
		b3_final = np.interp(z, red, b3)
		b4_final = np.interp(z, red, b4)
				
		return b1_final, b2_final, b3_final, b4_final

	elif bmodel == 3:
		for count,iz in enumerate(red):
			dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV'\
			'/case'+str(kcase)+'/coeff_3exp_0.0_z='+str(iz)+'.txt')
			with open(dat_file_path,'r') as f: 
				for i, line in enumerate(f):
					if i == mbin: 
						b1[count] = float(line.split()[0])
						b2[count] = float(line.split()[1])
						b3[count] = float(line.split()[2])
						b4[count] = float(line.split()[3])
					
		# interpolate on z
		b1_final = np.interp(z, red, b1)
		b2_final = np.interp(z, red, b2)
		b3_final = np.interp(z, red, b3)
		b4_final = np.interp(z, red, b4)
				
		return b1_final, b2_final, b3_final, b4_final
