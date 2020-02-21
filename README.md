# BE_HaPPy
Bias Emulator for Halo Power Spectrum is a software designed to facilitate future large scale surveys analysis by providing an accurate, easy to use and computionally inexpensive method to compute the halo bias in the presence of massive neutrinos.


1) You can download the folder or clone it from git url

2) To import everywhere from your laptop i used sys.path.append('/path/to/directory/') in example.py

The code was written with python2 but it should be also compatible with python3. The utilisation is simple, the user provides a linear Power Spectrum of his/her choice and the package will compute a new Power Spectrum according to the chosen configuration.


Dependencies
-------------
1. numpy
2. scipy
3. sys
4. os

Input linear Power Spectrum
---------------------------
By default the linear power spectrum provided is the one provided by Fast-PT equivalent by P_11. Instead the user can import the linear power spectrum of his/her choice. Some examples are provided in example.py

For consistency it is better to Pcc or Pcb

Bias models
---------------
1. Linear
2. Polynomial
3. Perturbation theory

Redshift models
---------------
1. Kaiser
2. Scoccimarro
3. TNS or eTNS (only available for the third bias model)


Other available options
-----
1. Real or Redshift space
2. The total neutrino mass
3. velocity dispersion for Finger of God
4. choice of mass bin or scale array
