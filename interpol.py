import numpy as np
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

def interpol_pt(pt_term, kbound, z, red, lred):
	
	size_k = 350 #size of the pt text file
	kpt = np.zeros((350))
	pt = np.zeros((350,lred))
	pt_temp = np.zeros((len(kbound),lred))
	#~ pt_final = np.zeros((len(kbound),len(np.atleast_1d(z))))
	pt_final = np.zeros((len(kbound)))


	for count,iz in enumerate(red):
		dat_file_path = os.path.join(dir_path, 'coefficients/0.0eV'\
		'/PT_coeff/'+pt_term+'_'+str(iz)+'.txt')

		with open(dat_file_path,'r') as f:  
			line = f.readline()
			for index_k in range(size_k):
				kpt[index_k] = float(line.split()[0])
				pt[index_k,count] = float(line.split()[1])
				line = f.readline()
		# interpolate on kbound
		pt_temp[:,count] = np.interp(kbound, kpt, pt[:,count]) 
		
		# interpolate on z
		for i in range(len(kbound)):
			#~ pt_final[i,:] = np.interp(z, red, pt_temp[i,:])
			pt_final[i] = np.interp(z, red, pt_temp[i,:])
			
	return pt_final
	
