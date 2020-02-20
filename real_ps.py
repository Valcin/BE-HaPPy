import numpy as np
import os
import math

dir_path = os.path.dirname(os.path.realpath(__file__))

def real_ps(mbin, bmodel, kbound, b1, b2, b3, b4, A, B, C, D, E, F, G, H, Pmod_dd, alpha):

	if bmodel == 1:
		b = b1*alpha
		Preal = b**2*Pmod_dd
		
	elif bmodel == 2:
		b = (b1 + b2*(kbound**2) + b3*(kbound**3) + b4*(kbound**4)) * alpha
		Preal = b**2*Pmod_dd
		
	elif bmodel == 3:
		# here b3 == bs and b4 == b3nl
		PhhDD = b1**2*Pmod_dd + b1*b2*A + 1/4.*b2**2*B + b1*b3*C + 1/2.*b2*b3*D + 1/4.*b3**2*E +\
		2*b1*b4*F 
		Preal = PhhDD * alpha**2
	
	return Preal
