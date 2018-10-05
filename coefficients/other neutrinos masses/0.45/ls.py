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
Mnu       = 0.03  #eV
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


	nv = 0.45



	#----------------------------------------------------------------
	#----------Tinker, Crocce param bias ----------------------------
	#----------------------------------------------------------------

	#~ #compute tinker stuff
	limM = [4.2e11,1e12,3e12,1e13, 5e15]
	loglim = [ 11.623, 12., 12.477, 13., 15.698]
	cname = 'Pk_cc_z='+str(z[j])+'00.txt'
	camb = np.loadtxt(cname)
	kcamb = camb[:,0]
	Pcamb = camb[:,1]
	
	massf = np. loadtxt('/home/david/codes/Paco/data2/0.15eV/hmf_z='+str(z[j])+'.txt')
	m_middle = massf[:,10]
	dm = massf[:,11]

	#### get tinker and crocce hmf
	bt=np.empty(len(m_middle),dtype=np.float64)
	for i in range(len(m_middle)):
			bt[i]=bias(kcamb,Pcamb,Omega_m,m_middle[i],'Tinker')
			
			
	dndM=MFL.Tinker_mass_function(kcamb,Pcamb,Omega_m,z[j],limM[0],limM[4],len(m_middle),Masses=m_middle)[1]
	with open('/home/david/codes/Paco/data2/other neutrinos masses/'+str(nv)+'/thmf_z='+str(z[j])+'.txt', 'w+') as fid_file:
	#~ with open('/home/david/codes/Paco/data2/other neutrinos masses/'+"%.3" % nv+'/thmf_z='+str(z[j])+'.txt', 'w+') as fid_file:
		for m in xrange(0, len(m_middle)):
			fid_file.write('%.8g\n' % (dndM[m]))
	fid_file.close()
	dndMbis=MFL.Crocce_mass_function(kcamb,Pcamb,Omega_m,z[j],limM[0],limM[4],len(m_middle),Masses=m_middle)[1]
	with open('/home/david/codes/Paco/data2/other neutrinos masses/'+str(nv)+'/chmf_z='+str(z[j])+'.txt', 'w+') as fid_file:
	#~ with open('/home/david/codes/Paco/data2/other neutrinos masses/'+"%.2" % nv+'/chmf_z='+str(z[j])+'.txt', 'w+') as fid_file:
		for m in xrange(0, len(m_middle)):
			fid_file.write('%.8g\n' % (dndMbis[m]))
	fid_file.close()
		
		
	dndM = np.loadtxt('/home/david/codes/Paco/data2/other neutrinos masses/'+str(nv)+'/thmf_z='+str(z[j])+'.txt')
	dndMbis = np.loadtxt('/home/david/codes/Paco/data2/other neutrinos masses/'+str(nv)+'/chmf_z='+str(z[j])+'.txt')
	bin1 = np.where(m_middle > 5e11 )[0]
	bin2 = np.where(m_middle > 1e12 )[0]
	bin3 = np.where(m_middle > 3e12 )[0]
	bin4 = np.where(m_middle > 1e13 )[0]	
	
	### CROCCE
	Bias_eff_t1=np.sum(dndMbis[bin1]*dm[bin1]*bt[bin1])/np.sum(dm[bin1]*dndMbis[bin1])
	Bias_eff_t2=np.sum(dndMbis[bin2]*dm[bin2]*bt[bin2])/np.sum(dm[bin2]*dndMbis[bin2])
	Bias_eff_t3=np.sum(dndMbis[bin3]*dm[bin3]*bt[bin3])/np.sum(dm[bin3]*dndMbis[bin3])
	Bias_eff_t4=np.sum(dndMbis[bin4]*dm[bin4]*bt[bin4])/np.sum(dm[bin4]*dndMbis[bin4])
	### TINKER
	bias_eff_t1=np.sum(dndM[bin1]*dm[bin1]*bt[bin1])/np.sum(dm[bin1]*dndM[bin1])
	bias_eff_t2=np.sum(dndM[bin2]*dm[bin2]*bt[bin2])/np.sum(dm[bin2]*dndM[bin2])
	bias_eff_t3=np.sum(dndM[bin3]*dm[bin3]*bt[bin3])/np.sum(dm[bin3]*dndM[bin3])
	bias_eff_t4=np.sum(dndM[bin4]*dm[bin4]*bt[bin4])/np.sum(dm[bin4]*dndM[bin4])

	with open('with_Crocce_HMF/LS_z='+str(z[j])+'_.txt', 'w+') as fid_file:
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (Bias_eff_t1,Bias_eff_t2, Bias_eff_t3, Bias_eff_t4))
	fid_file.close()
	with open('with_Tinker_HMF/LS_z='+str(z[j])+'_.txt', 'w+') as fid_file:
		fid_file.write('%.8g %.8g %.8g %.8g\n' % (bias_eff_t1,bias_eff_t2, bias_eff_t3, bias_eff_t4))
	fid_file.close()
