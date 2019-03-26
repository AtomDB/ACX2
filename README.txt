The AtomDB CX model is a model of charge exchange during a collision between a recombining charged ion and a donor atom or ion. As the electron is transferred from the donor to the recombining ion, it forms the recombined ion in an excited state. As this recombined ion relaxes to its ground state, it releases many photons.

An original version of ACX was released in 2016, which used empirical formulae for CX emission of all ions of all the elements up to nickel. These formulae, crucially, did not include any velocity dependent effects, which are important for correctly calculating the n, l and S of the excited levels captured into.

We have now taken data from the Kronos database ([1]_, [2]_, [3]_), which covers many fully stripped and one electron ions, and included it here. This has created a much improved dataset, which correctly captures the energy dependence of the process for these ions. For other ions not in the Kronos, the data falls back on ACX1 behaviour.

Once Kronos or ACX1 have been used to calculated the correct capture cross sections for each n, l and/or S shell, the data is combined with the AtomDB database (www.atomdb.org) to calculate the cascade path to ground, and the subsequent emissivities and wavelengths.

Included in this package are 2 additional features:

.. [1] Mullen, P. D., et al. ApJS 224, 31 (2016)
.. [2] Mullen, P. D., et al. ApJ 844, 7 (2017)
.. [3] Cumbee, R. S., et al. ApJ 852, 7 (2018)

=======================
Installation
=======================
Standard python installation:
python setup.py install

=====
Usage
=====

There are several classes in the acx2.py file. For most people, these will be irrelevant, and should be ignored. The main ones depend on the use case.

-----
XSPEC
-----

This model only works with the python XSPEC interface, pyxspec_. 
For xspec users, the acx2_xspec module contains what you need. From a python3 shell, and with xspec loaded ("import xspec"), run "import acx2_xspec" and all of the data will load. By default, these files contain the H and He donor files, as does the ACX module.

Three different models are loaded:

- acx2 : Emission from CX with the 14 main elements. Abundance is tied between all elements (so there is only 1 abundance keyword). Analogous to the apec model.
- vacx2 : Emission from CX with the 14 main elements. Abundance is free to vary between all the elements (though it starts frozen). Analagous to the vapec model.
- vvacx2 : Emission from 27 elements, H through Ni excluding Co. Abundance is free to vary between all the elements.

Note that in the acx and vacx cases, unlike in the apec and vapec models, the effective abundance of the minor recombining elements is 0, not solar. This speeds up calculation time and does not significantly effect the resulting emission.

++++++++++++++++
Model parameters
++++++++++++++++

+--------------+-----------------------------------------------------------------------------------+
| Parameter    | Definition                                                                        |
+==============+===================================================================================+
| temperature  | Plasma temperature (keV). Used for recombining particle ion fraction              |
+--------------+-----------------------------------------------------------------------------------+
| collnpar     | Collsion parameter (kev/u,km/s). Reduced energy or velocity of collision          |
+--------------+-----------------------------------------------------------------------------------+
| collntype    | Sets meaning of collnpar:                                                         |
+--------------+-----------------------------------------------------------------------------------+
|              |  1 - reduced energy (kev/u)                                                       |
+--------------+-----------------------------------------------------------------------------------+
|              | 2 - center of mass velocity (km/s)                                                |
+--------------+-----------------------------------------------------------------------------------+
|	       | 3 - donor ion velocity (km/s)                                                     |
+--------------+-----------------------------------------------------------------------------------+
|	       | 4 - recombining ion velocity (km/s)                                               |
+--------------+-----------------------------------------------------------------------------------+
| acxmodel     | ACX model to fall back on, from 1 to 8.                                           |
+--------------+-----------------------------------------------------------------------------------+
| recombtype   | single recombination (1) or all the way to neutral (2)                            |
+--------------+-----------------------------------------------------------------------------------+
| Hefrac       | Number fraction of donor which is He (remainder is H).                            |
+--------------+-----------------------------------------------------------------------------------+
| abund        | recombining elemental abundances. (given by individual element in vacx and vvacx) |
+--------------+-----------------------------------------------------------------------------------+

===============
Version History
===============
0.1.0
March 15th 2019
Initial release



.. _pyxspec: https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html
