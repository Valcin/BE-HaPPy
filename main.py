from interpol import interpol_pt
from coeff import bcoeff
from rescaling import rescaling
from real_ps import real_ps
from power_spec import red_ps

def ps_calc(coord,kcase, Mnu, mbin, rsd, bmodel, kbound, z, fz, Dz, fog, sigma_v = None, Plin = None):
	print('you chose: ')
	if coord == 0:
		print('- real space')
	elif coord == 1:
		print('- redshift space')
	print('- Total neutrino mass = ' + str(Mnu)+'eV')
	print('- the mass bin M' + str(mbin +1))
	if rsd == 1:
			print( '- the Kaiser model')
	elif rsd == 2:
			print('- the Scoccimaro model')
	elif rsd == 3:
			print( '- the TNS model')
	if bmodel == 1:
		print('- the linear bias')
	elif bmodel == 2:
		print('- the polynomial bias')
	elif bmodel == 3:
		print('- the perturbation theory bias')
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
	
	if Plin == None:
		Plin = Pmod_dd
	
	####################################################################
	####################################################################
	### load bias coefficients
	
	b1, b2, b3, b4 = bcoeff(mbin, bmodel, z, red, lred, kcase)
	
	####################################################################
	####################################################################
	### compute the neutrino rescaling coefficient

	alpha  = rescaling(z, mbin, Mnu)
	
	####################################################################
	####################################################################
	### compute the redshift power spectrum
	
	if coord == 0:
		Power = real_ps(mbin, bmodel, kbound, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Plin, alpha)
	elif coord == 1:
		Power = red_ps(mbin, bmodel, kbound, z, fz, Dz, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Plin,
	Pmod_dt, Pmod_tt, alpha, fog, rsd, red, kcase)

	return Power
