The AtomDB CX model is a model of charge exchange during a collision between a recombining charged ion and a donor atom or ion. The electron is transferred from the donor to the recombining ion, forming a recombined ion often in an excited state. As this recombined ion relaxes to its ground state, it releases a cascade of photons with relative intensities charactersitic of CX recombination.

An original version of ACX was released in 2016, which used empirical formulae for CX emission of all ions of all the elements up to nickel. These formulae, crucially, did not include any velocity dependent effects, which are important for correctly calculating the n, l and S of the excited levels captured into. In addition, spectral information was hardwired and difficult to update, resulting in updates to AtomDB not often being reflected in the following charge exchange spectra.

We have now taken CX cross section data from the Kronos database ([1]_, [2]_, [3]_), which covers many fully stripped and one electron recombining ions, and included it here. This has created a much improved dataset, which correctly captures the energy dependence of the process for these ions. For other ions not in the Kronos database, the model falls back on ACX1 behaviour.

Once Kronos or ACX1 have been used to calculate the correct capture cross sections for each n, l and/or S shell, the data is combined with the AtomDB database (www.atomdb.org) to calculate the cascade path to ground, and the subsequent emissivities and wavelengths. For ions with capture in to highly excited levels which AtomDB doesn't contain, AUTOSTRUCTURE calculations are preformed to get energy levels, wavelength and A-values for these transtitions. The result is a set of 3 files for each donor ion. The sigma file contains the cross section information for each ion. The line and cont[inuum] files contain the line emission and continuum emission from each shell capture in to, with a resolution appropriate for the model in question. For example, for ions with *nlS* resolved Kronos data, a spectrum is produced for each *n*, *l* and *S* capture and subsequent cascade. Thus there can be numerous spectra for each ion - there are 239 entries for Cl\ :sup:`7+`\ reflecting each *n*, *l* and *S* which Kronos contains cross sections for. For ACX-level data, the spectra are calculated for each relevant *n* and the four *l* distributions, as outlined in the ACX documentation (even, statistical, Landau-Zener and separable).

For a given interaction velocity or energy, the model uses the center of mass energy to obtain the cross section for capture into each shell from Kronos. For each shell where the cross section is greater than zero, a spectrum is calculated from the *line* and *cont* files. These are then multiplied by the appropriate cross section and summed to give the spectrum for CX of a particular ion.

.. [1] Mullen, P. D., et al. ApJS 224, 31 (2016)
.. [2] Mullen, P. D., et al. ApJ 844, 7 (2017)
.. [3] Cumbee, R. S., et al. ApJ 852, 7 (2018)

=======================
Installation
=======================
Standard python installation

.. code-block:: bash

   python setup.py install

.. note::
   ACX2 is a python 3 only module. Depending on your system's setup, you may need to substitute ``python3`` for all references to ``python``.

There are several useful flags that can be provided to this call, depending on your system:

  - ``--user`` flag causes installation in the user's home directory (useful if you lack root priviliges)

  - ``develop`` instead of ``install`` will install links to the current directory. This is useful if you want to edit/debug/develop the files further.

=====
Usage
=====


Each model in ACX2 can have an arbitrary set of donors. By default for the XSPEC model these are neutral H and He, but others may be selected. Additional input files will be required for these - please contact the project via the AtomDB or GitHub pages to make or discuss requests.

----------
Data Files
----------

Each model requires a set of data files to be installed with it. As these files are large they cannot be exported through GitHub. These files will be downloaded automatically to your $ATOMDB folder on first opening the package.

The files for each donor are:
  - ``sigma`` files: the cross sections for capture into each n, l and S (depending on the ion) from the Kronos database
  - ``line`` files: the line emission for capture into each n, l and S (depending on the ion) or each n and ACX1 l distribution for ions with no Kronos data
  - ``cont`` files: same as line files, but including continuum emission. True continuum in CX is entirely 2-photon emission from H-, He- and Be-like ions.

.. note::
  The emissivity data files have thousands of HDUs as currently assembled. Although these files read quickly in python, when opening in some programs (e.g. ``fv``) the load times can be upwards of 10 minutes. **As of January 2024 this (version 2.1) this has been resolved.

-------
Classes
-------
The ``acx2.py`` file contains a range of classes which can be used to model different aspects of the charge exchange. The basic principal is that the fits files contain the emissivity for each ion, broken down to reflect the way that Kronos handles the data.

+------------------+-------------------------+-------------------------------------------------+
|Kronos Resolution | Typical recombining ion | ACX2 handling                                   |
+==================+=========================+=================================================+
|n, l, S resolved  | hydrogenic              | Capture into each n, l, S                       |
|                  | bare C,N,O,Ne           |                                                 |
+------------------+-------------------------+-------------------------------------------------+
|n resolved        | bare                    | Capture into each n, ACX for l distribution     |
+------------------+-------------------------+-------------------------------------------------+
|not included      | all others              | Capture into 2 n shells, ACX for l distribution |
+------------------+-------------------------+-------------------------------------------------+

To handle this, the acx2 module contains 4 levels of classes:

- ``ACXModel`` : The overall ACX model. Can include multiple donor ACXDonorModel objects.

- ``ACXDonorModel`` : The ACX model for one donor. Contains spectra from each recombining ion.

- ``CXIonSpectrum`` : The spectrum for one recombining ion. Placeholder for ``CXIonSpectrum_ACX1``, ``CXIonSpectrum_N``, ``CXIonSpectrum_NLS`` classes, which handle the 3 cases is the table above. Contains CXShellSpectrum as required to get the data.

- ``CXShellSpectrum`` : The actual spectrum from a single each n, l, S shell (or n, ldist shell).


-----
XSPEC
-----

To use the model in XSPEC, one can ignore the class details above. Unfortunately, the code only works with the XSPEC python interface, pyxspec_ for now. Before loading the code, you will need to edit the ``acx2_xspec.py`` file to change the data file paths.

.. note::
  You will need to edit the acx2_xspec.py file:
  #. It may need to be moved into your path (depending on the data)
  #. The data file locations are hardcoded to $ATOMDB, you will need to update them to reflect where you have installed the line, continuum and cross section files if you have put them somewhere else.


To load the ACX2 model into XSPEC, acx2_xspec module contains what you need. From a python3 shell:

.. code-block:: python

  # import the xspec python module
  import xspec
  # import  acx2 wrapper
  import acx2_xspec

Once this is done, the data will load.


Five different models are loaded:

- acx2 : Emission from CX with the 14 main elements. Abundance is tied between all elements (so there is only 1 abundance keyword). Analogous to the apec model.
- vacx2 : Emission from CX with the 14 main elements. Abundance is free to vary between all the elements (though it starts frozen). Analagous to the vapec model.
- vvacx2 : Emission from 27 elements, H through Ni excluding Co. Abundance is free to vary between all the elements.
- acx2oneion : Emission from CX with one specific ion, specified with element and ion, where ion is the ion charge plus 1. 
- vacxnei : As with vacx2, but with the initial ionization balance set from a non-equilibrium plasma


.. note::
  Note that in the acx and vacx cases, unlike in the apec and vapec models, the effective abundance of the minor recombining elements is 0, not solar. This speeds up calculation time and does not significantly effect the resulting emission.

Once you have this, models can be used in pyxspec in the usual way, e.g.

.. code-block:: python

  m = xspec.Model('tbabs(pow+vacx2)')


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
|              | 1 - center of mass energy (kev/u)                                                 |
+--------------+-----------------------------------------------------------------------------------+
|              | 2 - center of mass velocity (km/s)                                                |
+--------------+-----------------------------------------------------------------------------------+
|              | 3 - donor ion velocity (km/s)                                                     |
+--------------+-----------------------------------------------------------------------------------+
|              | 4 - recombining ion velocity (km/s)                                               |
+--------------+-----------------------------------------------------------------------------------+
| acxmodel     | ACX model to fall back on, from 1 to 8.                                           |
+--------------+-----------------------------------------------------------------------------------+
| recombtype   | single recombination (1) or all the way to neutral (2)                            |
+--------------+-----------------------------------------------------------------------------------+
| Hefrac       | Number fraction of donor which is He (remainder is H).                            |
+--------------+-----------------------------------------------------------------------------------+
| abund        | recombining elemental abundances. (given by individual element in vacx and vvacx) |
+--------------+-----------------------------------------------------------------------------------+

.. note::
   The units for collision velocity in XSPEC are km/s, not cm/s as in the underlying ACX models. This is to keep the numbers closer to 1, which XSPEC likes.

++++++++++++++++++++
acxmodel definitions
++++++++++++++++++++
These are based on analytical formulae. For ions with nlS resolved
cross sections, these settings are ignored. For those with n resolved,
the l distribution is implemented, but the n distribution is ignored
(so models 1 and 5 give the same results). For those with no other
information, capture is into n shells defined by:

.. math::

  n' = q \sqrt{{{I_H}\over{I_p}}}\Big(1 + {{q-1}\over{\sqrt{2q}}}\Big)^{-1/2}
 
Where :math:`I_H` is the Rydberg constant, :math:`I_p` is the donor
ionization potential, and :math:`q` is the recombining ion charge.

For the acx1 models, the total cross section is
always fixed to :math:`3\times 10^{-15}` cm:sup:`2`. Capture is split
between the n shells, depending on the setting: if it is in one shell,
then this is the one closest to n' (rounding as normal). If it is
weighted, the split is (n'-n[round down]) into the n'[round up] shell,
and vice versa. So if n' is 6.25, then 1/4  and 3/4 of the emission goes into
the n=7 and 6 shells resepectively. 

+-------+-------------------+----------------------------------+
| Value | n distribution    | l distribution                   |
+=======+===================+==================================+
|  1    | one n shell       | even distribution by  l.         |
+-------+-------------------+----------------------------------+
|  2    | one n shell       | statistical distribution by l.   |
+-------+-------------------+----------------------------------+
|  3    | one n shell       | Landau-Zener distribution by  l. |
+-------+-------------------+----------------------------------+
|  4    | one n shell       | Separable distribution by l.     |
+-------+-------------------+----------------------------------+
|  5    | weighted 2 shells | even distribution by  l.         |
+-------+-------------------+----------------------------------+
|  6    | weighted 2 shells | statistical distribution by l.   |
+-------+-------------------+----------------------------------+
|  7    | weighted 2 shells | Landau-Zener distribution by l.  |
+-------+-------------------+----------------------------------+
|  8    | weighted 2 shells | Separable distribution by l.     |
+-------+-------------------+----------------------------------+


++++++++++++++++++++++++++++++++++++
Meaning of the l-shell distributions
++++++++++++++++++++++++++++++++++++

Note that all of these weighting schemes refer to the l distribution. For example, the
realtive weighting of total capture into levels with the 3s, 3p or 3d
configurations. Within all the levels sharing a configuration,
weighting is statistical (so 5 times more is captured into the :math:`1s3p ^3P_2` state than the :math:`1s3p ^{3}P_0)`.

+--------------+----------------------------------------------------------------------------------------------------------------+
| Value        | l distribution                                                                                                 |
+==============+================================================================================================================+
| Even         | Weighted evenly between l shells                                                                               |
+--------------+----------------------------------------------------------------------------------------------------------------+
| Statistical  | Weighted by the statistical weight of each level   .                                                           |
+--------------+----------------------------------------------------------------------------------------------------------------+
| Landau-Zener | Weighted by the function :math:`W(l)={{l(l+1)(2l+1)\times(n-1)! \times (n-2)!}\over{ (n+l)! \times (n-l-1)!}}` |
+--------------+----------------------------------------------------------------------------------------------------------------+
| Separable    | Weighted by the function :math:`W(l)={{(2l+1)}\over{Z}}\times \exp\Big[{{-l \times(l+1)}\over{z}}\Big]`        |
+--------------+----------------------------------------------------------------------------------------------------------------+


++++++++++++++++++++++++++
Normalization of the model
++++++++++++++++++++++++++

Ths model deals with two emissivities, which can get confusing. The photon emissivity of a line :math:`i\rightarrow j` is defined as:

.. math::
 \epsilon_{ij}=N_iA_{ij}

That is, the number of emitted photons :math:`cm^{-3} s^{-1}` is the number density :math:`N_i` of ions in state :math:`i`, times the spontaneous transition probability :math:`A_{ij}`.

We can ease the calculation of :math:`N_i` by separating out the calculation:

.. math::
  N_i = \frac{N_i}{N_{z1}}\frac{N_{z1}}{N_Z}\frac{N_Z}{N^r_H}N^r_H

where :math:`N_{z1}` is the ion abundance and :math:`N_{Z}` is the element abundance. The :math:`N_{z1}/N_{Z}` term is set by the temperature parameter, which is used to set the ion fraction. The :math:`N_{Z}/N^r_{H}`  is set for the recombining plasma by the abundance parameter, relative to the solar values of Anders and Grevesse 1989.

The ACX2 model solves the :math:`N_i/N_{z1}` problem by setting up a radiative matrix, with levels populated by CX and then radiative decay to the ground state forming the rest of the matrix. The CX rate coefficient into level :math:`i` is given by

.. math::
  \alpha^{CX}_{i} (cm^{3}s^{-1}) = \langle v_{com} \sigma_{i}(E)\rangle

and the rate per recombining ion per second is

.. math::
  \alpha^{CX}_{i} (s^{-1}) = \langle v_{com} \sigma_{i}(E)\rangle\ N^d_{H}

We solve the radiative matrix to obtain :math:`N_i/N_{z1}`, without the donor densities included (as they are multipliers on all the diagonal matrix elements we are effectively just moving them outside the matrix). This leaves us with:

.. math::
  \epsilon_{ij} = \frac{N_i}{N_{z1}}\frac{N_{z1}}{N_Z}\frac{N_Z}{N_H}N^r_HA_{ij}  N^d_{H}

ACX2 calculates the photon emissivity coefficient, :math:`\varepsilon_{ij}=\frac{N_i}{N_{z1}}A_{ij}`, and multiplies in the elemental and ion abundances based on the abundance and temperatures specified. This leaves:

.. math::
  \epsilon_{ij} = \varepsilon_{ij} \frac{N_{z1}}{N_Z}\frac{N_Z}{N_H} N^r_H N^d_{H}

  \epsilon_{ij} = \left(\mathrm{ACX2output}\right) N^r_H N^d_{H}


To convert this to a flux from a source to our instrument we integrate over the emitting volume and account for radiation over :math:`4\pi`. We also, at this point, repeat the process for the He donor and add the results, accounting for the different donor ion fractions.

.. math::

  \Gamma_{ij} (cm^{-2}s^{-1}) =  \frac{\int (\mathrm{ACX2output}) N^r_H N^d_{(H+He)} dV} {4 \pi D^2}


++++++++++++++++++++++++++++++++
XSPEC Normalization of the model
++++++++++++++++++++++++++++++++

The geometric norm represents the geometric parts of the flux calculation with a single number:

.. math::

  \text{norm}_{\text{geom}} (cm^{-5}) =  \frac{\int N^r_H N^d_{(H+He)} dV} {4 \pi D^2}


The XSPEC normalization is adjusted from the above in a few ways to make fitting more reliable. First, there is a factor of :math:`10^{10}` applied to bring the value closer to 1, which makes XSPEC fitting more reliable.

Secondly, for versions :math:`\geq 1.1.0`, the normalization is divided by the center of mass velocity. This has been implemented to compensate for the increase in flux with an increase in velocity (since :math:`\varepsilon_{ij} \propto \langle v_{com} \sigma_{i}(E)\rangle`), which resulted in the norm being anticorrelated with the collnpar. As this value would be different for every ion, the correction factor is based on a carbon-12 recombining ion and a hydrogen donor.

To recover the true emissivity of the plasma given an XSPEC fit result:

#.  If using version :math:`\geq 1.1.0`, calculate the correction factor, cf:

    #. If collntype == 1: cf = numpy.sqrt(4786031.3*collnpar/25.)
    #. If collntype == 2: cf = 1.0 * collnpar
    #. If collntype == 3: cf = 1.0 * collnpar/(1.0+12.0)  = collnpar/13
    #. If collntype == 4: cf = 12.0 * collnpar/(1.0+12.0)  = collnpar*12/13

#.  Else, cf = 1
#.  :math:`\text{norm}_{\text{geom}} = \text{norm}_{\text{XSPEC}} * \text{cf} * 10^{-10}`



===============
Version History
===============
1.0.0
March 15th 2019
Initial release

1.0.1
October 25th 2019
Fixed error in vacx2 XSPEC interface, which specified but did not implement fluorine
leading to an off-by-one error for all higher-Z elements

1.0.2
February 27th 2020
Error in velocity unit conversion corrected, thanks to Gabrielle Betancourt-Martinez for reporting the bug. This will not have affected fits performed through XSPEC

1.0.3
July 9th 2020
Updated code for compatibility with changes in the PyAtomDB interface

1.1.0
November 16th 2022
Major changes to the normalization. It now has the center of mass velocity of carbon-12 divided out of it.
This removed the velocity-normalization correlation which was otherwise present.

Added redshift to parameters.

Converted XSPEC interface collntype, acxmodel and recombtype into integer switches

1.1.1
December 1st 2022
Added extra option, 'calc_line_emissivity', which returns the emissivity of a specific transition due to CX. This can also be accessed during XSPEC sessions.
Put more examples in the new "examples" directory

1.1.2
December 2nd 2022
Bugfix to 'calc_line_emissivity', updates to examples.

1.1.3
January 18th 2023
Bugfix to 'acx2_xspec' interface: vacx and vvacx had incorrect parameter indexes

2.1.0
January 8th 2024
Version number incremented to sync data and code release file numbers
Major update to add velocity and thermal broadening. Data files rearranged to enable faster processing.
A big thank you to McKenna Blake for working on this project.

2.1.1
May 28th 2024
Bugfix: acx2_xspec has mis-indexed the redshift, leading to incorrect energies. Thank you to Shuinai Zhang for identifying this mistake.

2.1.2
June 16th 2024
Renamed redshift to Redshift in XSPEC wrapper for consistency with
other XSPEC models. Add in oneacx2 xspec model.

2.1.3
July 25th 2024
Added a third recombination option, where each ions' recombination is normalized by the total cross section.
Fixed bug with calc_line_emissivity where normalization was incorrect. Thank you to Yuki Amano for spotting this bug.

2.1.4
August 8th 2024
Replaced deprecated scipy.integrate.cumptrapz with cumulative_trapezoid

2.2.0
February 13th 2025
Added a new function to list the lines in a wavelength region and then get the quantum numbers. See examples/test_acx_linelist.py

2.3.0
April 7th 2025
Added a new model, vacxnei, for non-equilibrium plasma. See examples/test_acx_linelist.py

2.4.0
July 16th 2025
Reworked code to increase model speed, should be ~ 10 times faster (after initialization)

2.4.1
August 15th 2025
Automatic installation of required datafiles should now work.

.. _pyxspec: https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html
