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

Each model requires a set of data files to be installed with it. As these files are large they cannot be exported through GitHub, and they should instead be downloaded from the AtomDB CX webpage, www.atomdb.org/CX.

The files for each donor are:
  - ``sigma`` files: the cross sections for capture into each n, l and S (depending on the ion) from the Kronos database
  - ``line`` files: the line emission for capture into each n, l and S (depending on the ion) or each n and ACX1 l distribution for ions with no Kronos data
  - ``cont`` files: same as line files, but including continuum emission. True continuum in CX is entirely 2-photon emission from H-, He- and Be-like ions.

.. note::
  The emissivity data files have thousands of HDUs as currently assembled. Although these files read quickly in python, when opening in some programs (e.g. ``fv``) the load times can be upwards of 10 minutes. Rearranging these files to not cause this issue is a priority to fix.

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
  #. The data file locations are hardcoded, you will need to update them to reflect where you have installed the line, continuum and cross section files.


To load the ACX2 model into XSPEC, acx2_xspec module contains what you need. From a python3 shell:

.. code-block:: python

  # import the xspec python module
  import xspec
  # import  acx2 wrapper
  import acx2_xspec

Once this is done, the data will load.


Three different models are loaded:

- acx2 : Emission from CX with the 14 main elements. Abundance is tied between all elements (so there is only 1 abundance keyword). Analogous to the apec model.
- vacx2 : Emission from CX with the 14 main elements. Abundance is free to vary between all the elements (though it starts frozen). Analagous to the vapec model.
- vvacx2 : Emission from 27 elements, H through Ni excluding Co. Abundance is free to vary between all the elements.

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

===============
Version History
===============
0.1.0
March 15th 2019
Initial release



.. _pyxspec: https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html
