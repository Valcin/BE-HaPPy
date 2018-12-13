# BE_HaPPy
Bias Emulator for Halo Power Spectrum is a software designed to falicitate future large scale surveys analysis by providing an accurate, easy to use and computionally inexpensive method to compute the halo bias in the presence of massive neutrinos.



solutions to make the code available in every folder:
-----------------------------------------------------
1) install in the directory where you want to use the code
2) temporary solution sys.path.append('/path/to/directory/') in the python script where you want to import the code
3) export the directory where you installed BE-HaPPy in your .bashrc file

dependencies
-------------
need numpy
need scipy
explain how to define nuisance parameter
need mpk, mtk, halofit, 'z_max_pk': self.zmax

Redshift model
---------------
0. for kaiser
1. for scoccimarro
2. for tns
3. for etns (only available for exp is chosen)

fog
-----
veldisp range taken for hector gil marin Perturbation theory approach for the power spectrum: from dark matter in
real space to massive haloes in redshift space
euclid_pk.use_nuisance = ['P_shot','sigma_v']
data.parameters['sigma_v'] = [6, 0, 12,0.01,1,'nuisance']

my class
---------
expalin i put my version of class
