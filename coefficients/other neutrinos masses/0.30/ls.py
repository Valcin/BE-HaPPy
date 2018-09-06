import numpy as np
import redshift_space_library as RSL
from readfof import FoF_catalog
import MAS_library as MASL
import Pk_library as PKL
import mass_function_library as MFL
import bias_library as BL
from time import time
from bias_library import halo_bias, bias

z = [0.0,0.5,1.0,2.0]
#~ z = [0.0,2.0]
#~ mu = 0.5
#~ kmax = 1
#~ mass_range = ['m1','m2','m3','m4']
#~ mass_range = ['m1', 'm2']
#~ mass_range = ['m1']
#~ axis = 0 #in redshift-space distortion axis

# neutrino parameters
hierarchy = 'degenerate' #'degenerate', 'normal', 'inverted'
Mnu       = 0.30  #eV
Nnu       = 0  #number of massive neutrinos
Neff      = 3.046

# cosmological parameters
h       = 0.6711
Omega_c = 0.2685 - Mnu/(93.14*h**2)
Omega_b = 0.049
Omega_l = 0.6825
Omega_k = 0.0
Omega_m = Omega_c + Omega_b
tau     = None

#~ # read snapshot properties
#~ head = readsnap.snapshot_header(snapshot_fname)
BoxSize = 1000.0 #Mpc/h                                         
#~ redshift = head.redshift
#~ Hubble = 100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)#km/s/(Mpc/h)
#~ h = head.hubble

start = time()

for j in xrange(0,len(z)):
########################################################################
########################################################################
	####################################################################
	##### scale factor 

	red = ['0.0','0.5','1.0','2.0']
	ind = red.index(str(z[j]))
	f = [0.524,0.759,0.875,0.958]
	Dz = [ 1.,0.77,0.61,0.42]
	print 'For redshift z = ' + str(z[j])
	
	Omeg_m_z = Omega_m * (1 + z[j])**3 / (Omega_m * (1 + z[j])**3 + Omega_l)
	
	########################################################################
	############# 	0.0 eV Masseless neutrino 


	nv = 0.30


	#----------------------------------------------------------------
	#----------Tinker, Crocce param bias ----------------------------
	#----------------------------------------------------------------

	#~ #compute tinker stuff
	limM = [5e11,1e12,3e12,1e13, 3.2e15]
	limM = [4.5e11,1e12,3e12,1e13, 3.5e15]
	loglim = [ 11.653, 12., 12.477, 13., 15.544]
	cname = 'Pk_cb_z='+str(z[j])+'00.txt'


	#### get tinker and crocce hmf
	Tb1 = halo_bias('Tinker', z[j], limM[0],limM[4], cname,Omega_c, Omega_b, do_DM=True )
	Tb2 = halo_bias('Tinker', z[j], limM[1],limM[4], cname,Omega_c, Omega_b, do_DM=True )
	Tb3 = halo_bias('Tinker', z[j], limM[2],limM[4], cname,Omega_c, Omega_b, do_DM=True )
	Tb4 = halo_bias('Tinker', z[j], limM[3],limM[4], cname,Omega_c, Omega_b, do_DM=True )

	print Tb1, Tb2, Tb3, Tb4
	with open('LS_z='+str(z[j])+'_.txt', 'w+') as fid_file:
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (Tb1,Tb2, Tb3, Tb4))
	fid_file.close()
	
