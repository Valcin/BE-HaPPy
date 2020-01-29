import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

def rescaling(z, mbin, Mnu):
	
	dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV/large_scale/'\
	'LS_z='+str(z[0])+'_.txt')
	with open(dat_file_path,'r') as f: 
			line = f.readline()   
			bcc_LS000 = float(line.split()[mbin])
			
	
	nu_masses = [0.03, 0.06, 0.1, 0.13, 0.15, 0.3, 0.45, 0.6]
	bcc_massive = np.zeros((len(nu_masses)))
	for count, mn in enumerate(nu_masses):
		dat_file_path2 = os.path.join(dir_path, 'coefficients/other neutrinos masses/'+\
		str(mn)+'eV/LS_z='+str(z[0])+'_.txt')
		with open(dat_file_path2,'r') as f2: 
			line2 = f2.readline()   
			bcc_massive[count] = float(line2.split()[mbin])
			
	# interpolate on z
	bcc_final = np.interp(Mnu, nu_masses, bcc_massive)/bcc_LS000

	return bcc_final
