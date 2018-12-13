#### 1) copy the file "phalo.py" in your likelihood folder
#### 2) copy the file "bcc_coeff.dat" in the montepython data directory
#### 3) install python package halomod
#### 4) declare mTk in the need_cosmo_arguments of the __init__ file of your likelihood (below)
#### e.g. self.need_cosmo_arguments(data, {'output': 'mPk mTk'})
#### 4) request a non linear spectrum in the need_cosmo_arguments of the __init__ file of your likelihood
#### 5) declare a zmax for the computation of the transfer functions in the data file of your likelihood
#### e.g. self.need_cosmo_arguments(data, {'non linear': 'halofit'})
#### 6) you can import the halo power spectrum and the scale
#### e.g. from montepython.phalo import Phalo 
#### e.g. k, P_halo = Phalo(self,cosmo, data)
#### 6b) if there are "sub self" you can just modify the self in the function
#### e.g. for WiggleZ :  k, P_halo = Phalo(self.wigglez_a,cosmo, data)




1) Pk matter-real-space
2) Multipoles matter redshift-space
3) 2D Pk in redshift-space
4) save matter density field in real-space
5) save matter density field in redshift-space along 1 axis
6) save velocity field*
7) The halo catalogues can be saved without problems


Main File : phalo.py 
-----------------------
#### check if the transfer functions were computed (find a way to remove)

#### check if the non linear power spectrum has been requested (find a way to remove)

#### check if kmax is defined in the data file of your likelihood (find a way to remove)

#### For now the simulations available are for Mv=0,0.06,0.10,0.15 
#### Check if the total neutrino mass corresponds to one of the available ones
#### get_current_derived_parameters returns a dict so must be converted

#### import the requested redshift(s)

#### Store the selected redshifts in a array and deduce its length for the loops
#### array manipulation because len() and size only work for znumber >1

#### Store the redshifts where bcc fit  and bcc Ls are available in arrays
#### get the coefficients from the dat files in the data directory and do a linear interpolation on requested redshifts
#### get the large scale amplitude of bcc at different z and for different neutrino masses

#### Since get_transfer only gives transfer functions for specific k
#### and the cosmo.pk computes pk on given k, we must extract the k values from the get_transfer list 
#### so they coincide. 
#### For this "cosmo" they are 121 k values but we need to remove some terms of the arrays 
#### because they are out of bounds in in classy.Class.pk. CHECK WHY ???? (cf. euclid_lensing __init__ file)

#### get the index of the kmax value in the k.dot array and remove the elements out of range 
#### Since Raccanelli et al. 2017 are using a purely phenomenological fit over a reduced range of scales
#### if no kmax is defined the user is asked to define a kmax value in accordance with the 3 k-scale models in the article. 
#### This part of the code creates an array the values of kmax for each possible 
#### redshift 0,1,2 or 3.(Appox. values from Fig.4) and then select the k model chosen in the data file

#### select the k mode according to the ones in Raccanelli et al. 2017

#### create a scale array limited by kmin and kmax

#### import Omega_b and Omega_cdm from class. Remember to add Omega_cdm in classy and recompile after

#### get the value of h to rescale the power spectrum and wavenumber

#### define the CDM + baryons transfer function from Raccanelli et al 2017.

#### compute the bias using the fit formula from Raccanelli et al 2017.

#### get the non linear power spectrum from class and rescale it in (Mpc/h)**3

#### get the linear power spectrum from class

##### check the total neutrino mass and get the denom accordingly

#### rescale bcc_fit so the large scale amplitude matches the one from Paco (won't be needed after the correction)

#### compute the total matter bias bmm w.r.t bcc using formula 5 in Raccanelli et al.

#### Compute the halo Power spectrum in real space

#### Define the linear growth factor and growth rate (growth factor f in class)

#### Compute the matter Power spectrum monopole and quadrupole in redshift space + linear kaiser effect + (F.o.G effect TO ADD)

#### Compute the halo Power spectrum monopole and quadrupole in redshift space + linear kaiser  +(noise spectrum TO ADD)

#### Plot different parameters Pk, bcc_fit to check settings and write pk in a .dat file


RSD.py
----------------------

#### Query the cosmo module for the Hubble rate (in 1/Mpc) and convert it to km/s/Mpc
#### compute a velocity dispersion to simulate F.o.G effects

#### compute the redshift space distortion using the analytic expression of Linder F(k,mu,z)



























