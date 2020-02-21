# BE_HaPPy
Bias Emulator for Halo Power Spectrum is a software designed to falicitate future large scale surveys analysis by providing an accurate, easy to use and computionally inexpensive method to compute the halo bias in the presence of massive neutrinos.


1) You can download the folder or clone it from git url
2) To import everywhere from your laptop i used sys.path.append('/path/to/directory/') in example.py

The code was written with python2 but it should be also compatible with python3

dependencies
-------------
need numpy
need scipy
need sys
need os

Redshift models
---------------
1. for kaiser
2. for scoccimarro
3. for tns
4. for etns (only available for exp is chosen)

fog
-----
veldisp range taken for hector gil marin Perturbation theory approach for the power spectrum: from dark matter in
real space to massive haloes in redshift space
euclid_pk.use_nuisance = ['P_shot','sigma_v']
data.parameters['sigma_v'] = [6, 0, 12,0.01,1,'nuisance']
