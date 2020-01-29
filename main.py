from interpol import interpol_pt
from coeff import bcoeff
from rescaling import rescaling

def ps_calc(kcase, Mnu, mbin, rsd, bmodel, kbound, z):
	print 'Total neutrino mass is ' + str(Mnu)+'eV'
	print 'you chose the mass bin M' + str(mbin +1)
	if rsd == 1:
			print 'you chose the Kaiser model'
	elif rsd == 2:
			print 'you chose the Scoccimaro model'
	elif rsd == 3:
			print 'you chose the TNS model'
	if bmodel == 1:
		print 'you chose the linear bias'
	elif bmodel == 2:
		print 'you chose the polynomial bias'
	elif bmodel == 3:
		print 'you chose the perturbation theory bias'
	print ('')
	
	
	# store the calibrated redshift
	red = [0.0,0.5,1.0,2.0] 
	lred = len(red) 
	ind = red.index(z)
	
	####################################################################
	####################################################################
	### load perturbation terms 

	pt_terms = ['Pmod_dd', 'Pmod_dt', 'Pmod_tt','A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

	Pmod_dd = interpol_pt(pt_terms[0], kbound, z, red, lred)
	Pmod_dt = interpol_pt(pt_terms[1], kbound, z, red, lred)
	Pmod_tt = interpol_pt(pt_terms[2], kbound, z, red, lred)
	A = interpol_pt(pt_terms[3], kbound, z, red, lred)
	B = interpol_pt(pt_terms[4], kbound, z, red, lred)
	C = interpol_pt(pt_terms[5], kbound, z, red, lred)
	D = interpol_pt(pt_terms[6], kbound, z, red, lred)
	E = interpol_pt(pt_terms[7], kbound, z, red, lred)
	F = interpol_pt(pt_terms[8], kbound, z, red, lred)
	G = interpol_pt(pt_terms[9], kbound, z, red, lred)
	H = interpol_pt(pt_terms[10], kbound, z, red, lred)
	
	####################################################################
	####################################################################
	### load bias coefficients
	
	b1, b2, b3, b4 = bcoeff(mbin, bmodel, z, red, lred, kcase)
	
	####################################################################
	####################################################################
	### compute the neutrino rescaling coefficient

	alpha  = rescaling(z, mbin, Mnu)

	return
