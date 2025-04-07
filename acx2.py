"""
.. module:: acx2
   :platform: posix, macos
   :synopsis: Module for running the AtomDB Charge eXchange version 2

.. moduleauthor:: Adam Foster <afoster@cfa.harvard.edu>

"""

import pyatomdb, numpy, os, hashlib
import astropy.io.fits as pyfits

# this is the library for running the ACX2 model.

#some constants!

CENTEROFMASS=1
UNIVERSAL_CX_CROSSSECTION=3e-15 # in cm^2
DEBUG=False
SINGLE_RECOMBINATION = 1
FULL_RECOMBINATION = 2
FULL_NORMALIZED_RECOMBINATION = 3
__version__="2.1.4"

llistdtype = numpy.dtype({'names':['EPSILON','LAMBDA','ENERGY',\
                                   'ELEMENT','ION','ION_DRV',\
                                   'UPPERLEV','LOWERLEV'],\
                          'formats':[float, float, float, \
                                     int, int, int,\
                                     int, int]})


llistfulldtype = numpy.dtype({'names':['EPSILON','LAMBDA','ENERGY',\
                                       'ELEMENT','ION','ION_DRV',\
                                       'UPPERLEV','LOWERLEV',\
                                       'UPPERCFG','UPPERSYM',\
                                       'LOWERCFG','LOWERSYM',\
                                       'EINSTEIN_A', 'B_RATIO'],\
                              'formats':[float, float, float, \
                                         int, int, int,\
                                         int, int,\
                                         '|U40','|U40',\
                                         '|U15','|U15',\
                                         float, float]})

# This code is for initializing/auto-downloading the data files
import __init__

def loginterp(newx, x, y, offset = 1e-40):
  """
  Interpolation helper function. Interpolates linearly on a log-log grid
  If newx < x[0], return x[0].
  If newx > x[-1], extrapolate slope of last 2 points, unless y[-1] ==0
  in which case return 0.

  PARAMETERS
  ----------
  newx : float
    The new X parameters (energy, kev/amu)
  x : array(float)
    The x parameters (energy, kev/amu)
  y : array(float)
    The y parameters (cross section, cm2)
  offset : float (optional)
    An offset to apply before interpolation. Prevents  log(0) issues.

  RETURNS
  -------
  newy : float
    The interpolated cross section. Minimum of 0 in case of numerical issues.

  """

  import numpy


  if newx < x[0]:
    yout=y[0]

  else:

    if newx > x[-1]:
      if x[-1] > 0:
        yy = numpy.log(y+1e-40)
        xx = numpy.log(x+1e-40)
        newxx = numpy.log(newx+1e-40)

        yout = numpy.exp(((yy[-1]-yy[-2])/(xx[-1]-xx[-2])) *\
                (newxx-xx[-1])+yy[-1])-1e-40
      else:
        yout = 0.0
    else:
      yy = numpy.log(y+1e-40)
      xx = numpy.log(x+1e-40)
      newxx = numpy.log(newx+1e-40)

      yout = numpy.exp(numpy.interp(newxx, xx, yy))-1e-40
  return max([yout, 0.0])





class ACXModel():
  import pyatomdb, numpy, os, hashlib

  """
  This is the overarching model. Within it, it will load an
  ACXDonorModel for each donor atom/molecule. In theory, you should
  be able to load this in XSPEC and be done with it.

  PARAMETERS
  ----------
  None

  ATTRIBUTES
  ----------
  DonorList : list of ACXDonorModel
    List of donors, (e.g. H, He, H2...). Add to this list using add_donor
  ebins : ndarray(float)
    Energy bin edges for the resulting spectra
  ebins_checksum : md5sum
    The md5 hash of the energy bins. Stored to easily catch changes
  temperature : float
    The temperature in keV

  NOTES
  -----
  Once initialized, call "add_donor" to add donor elements.

  Provide energy bins using set_ebins, temperature with set_temperature

  Then call calc_spectrum to return total spectrum.
  """

  def __init__(self, elements=None):
    """
    Initialize the ACX Model

    elements: list(int)
      The elements to include as recombining ions (not the donors!)
      Defaults to H through Ni, excluding Co.
    """

    self.DonorList = []
    self.ebins = False
    self.ebins_checksum = False
    self.temperature = False
    self.ebins_set = False
    self.temperature_set = False
    self.collision_set = False
    self.recombtype = SINGLE_RECOMBINATION
    self.datacache={}


    if elements is None:
      self.elements = list(range(1,27))
      self.elements.append(28)
    else:
      self.elements = elements

    self.abund = {}
    for Z in self.elements:
      self.abund[Z] = 1.0


    print('AtomDB ACX version 2 model. Based on CX cross section data from the Kronos Database')
    print('See references:')
    print('  Mullen, P. D., et al. ApJS 224, 31 (2016)')
    print('  Mullen, P. D., et al. ApJ 844, 7 (2017)')
    print('  Cumbee, R. S., et al. ApJ 852, 7 (2018)')



  def add_donor(self, donor, \
                donor_linefile="$ATOMDB/acx2_XDONORX_vXVERSIONX_line.fits", \
                donor_contfile="$ATOMDB/acx2_XDONORX_vXVERSIONX_cont.fits", \
                donor_crosssectionfile="$ATOMDB/acx2_XDONORX_vXVERSIONX_sigma.fits", \
                abundset='AG89',\
                elements=None):
    """
    Add a donor to the ACX model (e.g. H$^{0+}$, He$^{0+}$)

    PARAMETERS
    ----------
    donor : str
      The donor symbol
    donor_linefile : str
      The file with the donor line emision
    donor_contfile : str
      The file with the donor continuum emision
    donor_crosssectionfile : str
      The cross section file for donor
    abundset : str
      The abundance set to use
    elements : list[int]
      The nuclear charge of the elements to recombine with.

    RETURNS
    -------
    None

    NOTES
    -----
    Sets parameters in self.DonorList
    """

    # Load default data values

    curversion = open(os.path.expandvars('$ATOMDB/VERSION_ACX'),'r').read()[:-1]

    #get the new filenames donor

    donor_linefile_actual = os.path.expandvars(donor_linefile).replace('XDONORX', donor).replace('XVERSIONX', curversion)
    donor_contfile_actual = os.path.expandvars(donor_contfile).replace('XDONORX', donor).replace('XVERSIONX', curversion)
    donor_crosssectionfile_actual = os.path.expandvars(donor_crosssectionfile).replace('XDONORX', donor).replace('XVERSIONX', curversion)

    if elements is None:
      elements=self.elements

    self.DonorList.append(ACXDonorModel(donor, donor_linefile_actual,\
                                        donor_contfile_actual,\
                                        donor_crosssectionfile_actual,\
                                        abundset=abundset,\
                                        elements=elements))
    if self.ebins_set:
      self.DonorList[-1].set_ebins(self.ebins)

    if self.temperature_set:
      self.DonorList[-1].set_temperature(self.temperature)

    if self.collision_set:
      self.DonorList[-1].set_collision(self.colltype, self.collunits)


  def set_acxmodel(self, acxmodel):
    """
    Set the ACX spectrum type

    PARAMETERS
    ----------
    acxmodel : int
      The acxmodel (between 1 and 8)

    """
    self.acxmodel=acxmodel

    for donor in self.DonorList:
      donor.set_acxmodel(self.acxmodel)


  def set_recombtype(self, recombtype):
    """
    Set recombination type: 1 = single, 2 = full, 3 = full renormalized


    PARAMETERS
    ----------
    recombtype : int
      The recombination type (1,2 or 3)

    """
    self.recombtype = recombtype

    for donor in self.DonorList:
      donor.set_recombtype(self.recombtype)


  # now set some spectral goodies
  def set_ebins(self, ebins):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """

    self.ebins = ebins
    self.ebins_checksum = hashlib.md5(ebins).hexdigest()
    self.ebins_set = True
    for donor in self.DonorList:
      donor.set_ebins(ebins)



  def set_temperature(self, temperature, teunits='keV'):
    """
    Set the ionization temperature for the data

    PARAMETERS
    ----------
    temperature : float
      The electron temperature
    teunits : str (K or keV)
      The electron temperature units

    RETURNS
    -------
    None

    """
    if teunits.lower() == 'k':
      self.temperature = temperature*pyatomdb.const.KBOLTZ
    elif teunits.lower() == 'kev':
      self.temperature = temperature
    else:
      print("Error, unknown temperature units in set_temperature")
      return

    self.temperature_set=True

    for donor in self.DonorList:
      if donor.temperature != self.temperature:
        # need to redefine
        donor.set_temperature(self.temperature)



  def set_ionfrac(self, ionfrac):
    """
    Recalculate the ionization balance based on equilibrium electron temperature

    PARAMETERS
    ----------
    ionfrac: dict of arrays
      ionization fraction, e.g. ionfrac[8]=numpy.array([0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.4, 0.2, 0.0])

    RETURNS
    -------
    None

    """

    # set the numbersset_collisionparam
    self.ionfrac_from_temperature = False
    self.temperature = False
    self.ionfrac = ionfrac

    # calculate ionization balance
    for donor in self.DonorList:
      donor.set_ionfrac(ionfrac)


  def set_donorabund(self, donorlist, abundlist):
    """
    Set the abundance of each donor (in DonorList[i].donor)

    PARAMETERS
    ----------
    donorlist : string or iterable of string
      donor symbols (e.g. H, He, H2O). Case insensitive.
    abundlist : float or array of float
      abundance of each.

    RETURNS
    -------
    None
    """


    try:
      _ = (d for d in donorlist)
    except TypeError:
      donorlist = [donorlist]
      abundlist = [abundlist]

    for i in range(len(donorlist)):
      donorlist[i] = donorlist[i].lower()


    for donor in self.DonorList:
      try:

        i = donorlist.index(donor.donor.lower())
        donor.set_donorabund(abundlist[i])
      except ValueError:
        print("no match for %s"%(donor.donor.lower()))
        pass



  def set_collisiontype(self, colltype, collunits):
    """
    Set the collision interaction frame of reference and values

    PARAMETERS
    ----------
    colltype: int {1,2,3,4}
      1 - energy of center of mass
      2 - velocity of center of mass
      3 - velocity of donor ion
      4 - velocity of recombining ion
    collunits : string {'cm/s', 'km/s', 'kev/amu','kev/u', 'ev/amu','ev/u'}
      The units of collvalue (case insensitive)
    """
#    print("Collunits", collunits)
    self.collisionunits=collunits
    self.collisiontype=colltype
    self.collisionset=True
    for donor in self.DonorList:
      donor.set_collisiontype(colltype, collunits=collunits)

  def calc_spectrum(self, collvalue, Tbroaden=0, vbroaden=0, ebroaden=0):
    """
    Calculate the spectrum for all the donors, sum.

    PARAMETERS
    ----------
    collvalue : float
      The collision parameter (kev/u or cm/s, depending) to calculate the spectrum
    Tbroaden : float
      The thermal broadening temperature (in K)
    vbroaden : float
      The velocity broadening (in km/s)

    RETURNS
    -------
    ret : array(float)
      The spectrum, in ph cm^3 s-1 bin-1
    """
    retset=False

    if DEBUG:
      ret = []
      for donor in self.DonorList:
        ret.append(donor.calc_spectrum(collvalue, Tbroaden, vbroaden))

      return ret
    else:
      ret=False
      for donor in self.DonorList:

        if retset==True:
          ret += donor.calc_spectrum(collvalue, Tbroaden, vbroaden, ebroaden)*donor.donorAbund
        else:
          ret = donor.calc_spectrum(collvalue, Tbroaden, vbroaden, ebroaden)*donor.donorAbund
          retset=True

      return ret

  def calc_linelist(self, collvalue, specrange, specunit='A', byDonor=False, fulldetail=False, redshift=0.0):
    """
    Calculate the spectrum for all the donors, sum.

    PARAMETERS
    ----------
    collvalue : float
      The collision parameter (kev/u or cm/s, depending) to calculate the spectrum
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)
    byDonor : bool
      If True, separate the results by donor. Otherwise sum identical lines from different donors
    fulldetail : bool
      If True, calls get_full_line_info to fill out the transition information
      e.g. configurations, quantum numbers.
    redshift : float
      Apply this redshift to the spectrum

    RETURNS
    -------
    linelist : numpy.array
      Numpy array of lines and emissivity, custom dtype. If byDonor=True, a dict with donor names and these data
   """
    retset=False

    if DEBUG:
      ret = []
      for donor in self.DonorList:
        ret.append(donor.calc_spectrum(collvalue, Tbroaden, vbroaden))

      return ret
    else:
      if byDonor:
        ret={}
      else:
        ret = numpy.zeros(0, dtype=llistdtype)

      for donor in self.DonorList:

        tmp = donor.calc_linelist(collvalue, specrange, specunit=specunit, redshift=redshift)
        tmp['EPSILON'] *=donor.donorAbund
        if byDonor:
          ret[self.donor] = tmp # e.g. ret['h'] is the collision with H
        else:
          ret=numpy.append(ret,tmp)


    if not byDonor:
        # need to filter/sum duplicates
      ret = _merge_linelist_duplicates(ret)

    return ret


  def set_abund(self, abund, elements=None):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    abund : array(float)
      Abundances, relative to defulat

    RETURNS
    -------
    None
    """

    try:
      if elements==None:
        elements = self.elements
      else:
        elements, elissvec = pyatomdb.util.make_vec(elements)
    except ValueError:
      elements, elissvec = pyatomdb.util.make_vec(elements)

    abundvec, aisvec = pyatomdb.util.make_vec(abund)

    if len(abundvec) != len(elements):
      if len(abundvec)==1:
        abundvec = abundvec[0]*numpy.ones(len(elements))

      else:
        raise(ValueError,"Error: specified incompatible number of elements from set_abund")
        return

    for i in range(len(abundvec)):
      self.abund[elements[i]] = abundvec[i]

    for donor in self.DonorList:
      donor.set_abund(abundvec, elements=elements)


  def calc_line_emissivity(self,  Z, z1, up, lo, collvalue=None):
    """
    Calculate the spectrum for all the donors, sum.

    PARAMETERS
    ----------
    collvalue : float
      The collision parameter (kev/u or cm/s, depending) to calculate the spectrum
    Z : int
      element charge
    z1 : int
      recombining ion charge +1.
    up : int
      upper level of transition
    lo : int
      lower level of transition

    RETURNS
    -------
    ret : dict
      Contains:
    Lambda : float
      The wavelength in Angstroms of the line
    Epsilon :  float
      The emissivity of the line in photons cm^3 s-1
    up : int
      The upper level of the line
    lo : int
      The lower level of the line


    """
    ret = {'Lambda': 0.0,
           'Epsilon': 0.0,
           'up': up,
           'lo': lo}
#    if DEBUG:
      #ret = []
      #for donor in self.DonorList:
        #ret.append(donor.calc_spectrum(collvalue, Tbroaden, vbroaden))

      #return ret
    #else:
    for donor in self.DonorList:

      lineemiss=  donor.calc_line_emissivity(Z, z1, up, lo, collvalue=collvalue)
      ret['Lambda'] = lineemiss['Lambda']
      ret['Epsilon']+=lineemiss['Epsilon']*donor.donorAbund

    return ret

  def calc_ionfrac_equilibrium(self):
    """
    Recalculate the ionization balance based on equilibrium electron temperature

    PARAMETERS
    ----------
    None

    RETURNS
    -------
    None

    NOTES
    -----
    Uses self.temperature (in keV) to set self.ionfrac
    """

    for Z in self.elements:
#        if not Z in self.ionfrac.keys():
          self.ionfrac[Z] = pyatomdb.atomdb.apec.return_ionbal(Z, self.temperature,\
                                                                   teunit='kev', datacache=self.datacache)
    self.ionfrac_from_temperature = True


  def calc_ionfrac_nonequilibrium(self, kT, kT_init, tau):
    """
    Recalculate the ionization balance based on equilibrium electron temperature

    PARAMETERS
    ----------
    kT : float
      electron temperature in keV
    kT_init : float
      initial ion distribution is equilibrium at this electron temperature (keV)
    tau : float
      ionization age (ne * t, cm^-3 s)

    RETURNS
    -------
    None

    NOTES
    -----
    Uses self.temperature (in keV) to set self.ionfrac
    """

    ionfrac = {}
    for Z in self.elements:
#        if not Z in self.ionfrac.keys():
          ionfrac[Z] = pyatomdb.apec.return_ionbal(Z, kT,\
                                                  init_pop=kT_init,\
                                                  tau = tau,\
                                                  teunit='kev', \
                                                  datacache=self.datacache,\
                                                  fast=True)
    #self.ionfrac_from_temperature = True
    self.set_ionfrac(ionfrac)



def _merge_linelist_duplicates(llist, keys=['ELEMENT', 'ION','UPPERLEV','LOWERLEV'], sumover='EPSILON'):
    """
    Look through the linelists and search for duplicates. Duplicated lines
    are summed.

    Parameters
    ----------
    llist : array(linelist)
      The list of lines to check
    keys : list
      Keys which must be identical to merge

    Returns
    -------
    llist_out : array(linelist)
      The lines after filtering for duplicates.
    """
    # sort the data
    llist =numpy.sort(llist, order=keys)
    # merge
    keep = numpy.ones(len(llist), dtype=bool)

    # find each level where the next is for the same transition

    if len(keys)==4:
      j = numpy.where((llist[1:][keys[0]]==llist[:-1][keys[0]]) &\
                      (llist[1:][keys[1]]==llist[:-1][keys[1]]) &\
                      (llist[1:][keys[2]]==llist[:-1][keys[2]]) &\
                      (llist[1:][keys[3]]==llist[:-1][keys[3]]))[0]
    elif len(keys)==3:
      j = numpy.where((llist[1:][keys[0]]==llist[:-1][keys[0]]) &\
                      (llist[1:][keys[1]]==llist[:-1][keys[1]]) &\
                      (llist[1:][keys[2]]==llist[:-1][keys[2]]))[0]

    for jj in j:
      # add the emissivity to the second of the 2
      llist[sumover][jj+1] += llist[sumover][jj]
      keep[jj]=False

    # remove all the non-keepers
    llist_out = llist[keep]

    return llist_out

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#--------BEGIN ACXDonorModel class----------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------




class ACXDonorModel():
  """
  A model of the ACX Donor

  PARAMETERS
  ----------
  donor : string
    The donor symbol (e.g. "H")
  donor_linefile : string
    The file with the donor line emision
  donor_contfile : string
    The file with the donor continuum emision
  donor_crosssectionfile : string
    The cross section file for donor
  abundset : optional, {"AG89"}
    The abundance set. Only AG89 works currently
  elements : optional, array-like int
    The recombining element atomic numbers. Default is all.
  acxmodel : optional, int, default=8
    The acx model l, n distribution to fall back on.

  ATTRIBUTES
  ----------
  linedata : HDUList
    The data in the donor_linefile
  contdata : HDUList
    The data in the donor_contfile
  crosssectiondata : HDUList
    The data in the donor_crosssectionfile
  donorAbundance : float
    The donor abundance. Defaults to 1.0. Intended for correctly mixing
    donors, e.g. in 10%He, 90%H donor plasma, set to 0.1 and 0.9 respectively,
    though no normalization is  checked for.
  donormass : float
    mass of the donor in a.m.u.
  abund : dict of float
    abundance of each element, relative to abundset. e.g. abund[6] = 2.0 means 2x solar C
  ionfrac_from_temperature : bool
    was the ionfracance calculated from the temperature
  temperature : float
    the temperature in keV
  ionfrac : dict of arrays of float
    the ionization fraction, normalized to 1 for each element,
    e.g. ionfrac[2] = ndarray([0.01,0.1,0.89])


  """
  import pyatomdb, numpy, os, hashlib

  def __init__(self, donor, \
               donor_linefile=None, \
               donor_contfile=None, \
               donor_crosssectionfile=None, \
               abundset='AG89',\
               elements=None,\
               acxmodel = 8, \
               recombtype = SINGLE_RECOMBINATION,\
               collisiontype = 1):

    if ((donor_linefile is None) or (donor_contfile is None) or (donor_crosssectionfile is None)):
      curversion = open(os.path.expandvars('$ATOMDB/VERSION_ACX'),'r').read()[:-1]
    if donor_linefile is None:
      donor_linefile_actual = os.path.expandvars("$ATOMDB/acx2_XDONORX_vXVERSIONX_line.fits").replace('XDONORX', donor).replace('XVERSIONX', curversion)
    else:
      donor_linefile_actual = os.path.expandvars(donor_linefile)

    if donor_contfile is None:
      donor_contfile_actual = os.path.expandvars("$ATOMDB/acx2_XDONORX_vXVERSIONX_cont.fits").replace('XDONORX', donor).replace('XVERSIONX', curversion)
    else:
      donor_contfile_actual = os.path.expandvars(donor_contfile)

    if donor_crosssectionfile is None:
      donor_crosssectionfile_actual = os.path.expandvars("$ATOMDB/acx2_XDONORX_vXVERSIONX_sigma.fits").replace('XDONORX', donor).replace('XVERSIONX', curversion)
    else:
      donor_crosssectionfile_actual = os.path.expandvars(donor_crosssectionfile)


    self.donor = donor.lower() # store in lower case
    self.donor_linefile = donor_linefile_actual
    self.donor_contfile = donor_contfile_actual
    self.donor_crosssectionfile = donor_crosssectionfile_actual


    try:
      self.linedata = pyfits.open(self.donor_linefile)
    except:
      print("Cannot open line data file %s"%(self.donor_linefile))
      raise

    try:
      self.contdata = pyfits.open(self.donor_contfile)
    except:
      print("Cannot open continuum data file %s"%(self.donor_contfile))
      raise

    try:
      self.crosssectiondata = pyfits.open(self.donor_crosssectionfile)
      self.donormass = self.crosssectiondata[1].header['DonMass']

    except:
      print("Cannot open cross section data file %s"%(self.donor_crosssectionfile))
      raise



    # create a structure for the spectral data
    self.spectra = {}

    # set the default abundance of the donor to 1.0
    self.set_donorabund(1.0)


    # create an abundance vecotr

    self.abundset=abundset
    self.default_abundset=abundset


    # if elements are specified, use them. Otherwise, use Z=1-30
    if elements is None:
      self.elements = list(range(1,27))
      self.elements.append(28)
    else:
      self.elements = elements


    # set the abundances:
    #   (1) the initial vector is whatever set AtomDB was calculated on,
    #       and is therefore 1.0
#    input('DON2.1')

    self.abundsetvector = {}
    for Z in self.elements:
      self.abundsetvector[Z] = 1.0
#    input('DON2.2')

    #   (2) but if another vector was already specified, use this instead
    if pyatomdb.util.keyword_check(abundset):
      self.set_abundset(abundset)
#    input('DON2.3')

    self.abund = {}
    for Z in self.elements:
      self.abund[Z]=1.0
#    input('DON2.4')

    # set up the default structure for the ionization fraction
    self.ionfrac = {}
    for Z in self.elements:
      self.ionfrac[Z]=numpy.zeros(Z+1, dtype=float)
      # default to everything fully stripped
      self.ionfrac[Z][Z] = 1.0
#    input('DON3')

    self.ionfrac_from_temperature = False
    self.temperature = False

    self.acxmodel= acxmodel

    self.set_collisiontype(collisiontype)

    # temporary things
    self.datacache={}

    self.ebins = 0.0
    self.ebins_checksum = ''

    self.recombtype = recombtype


  def set_donorabund(self, abund):
    """
    Set the donor abundance

    PARAMETERS
    ----------
    abund : float
      The abundance of the donor.

    """
    self.donorAbund = abund


  def set_abundset(self, abundstring):
    """
    Set the abundance set.

    Parameters
    ----------
    abundstring : string
      The abundance string (e.g. "AG89", "uniform"). Case insensitive.
      See atomdb.get_abundance for list of possible abundances

    Returns
    -------
    none
      updates self.abundset and self.abundsetvector.
    """
    # read in the abundance the raw data was calculated on
    old = pyatomdb.atomdb.get_abundance(abundset=self.default_abundset)


    # read in the new abundance
    new = pyatomdb.atomdb.get_abundance(abundset=abundstring)


    # divide the 2, store the replacement ratio to self.abundsetvector
    for Z in list(self.abundsetvector.keys()):
      self.abundsetvector[Z]=new[Z]/old[Z]


    # update the current abundance string to represent your input
    self.abundset=abundstring



  def find_crosssection_type(self, Z, z1):
    """
    Find the cross section type, to assign correct CXIonSpectrum type

    PARAMETERS
    ----------
    Z : int
      element charge
    z1 : int
      recombining ion charge +1.

    RETURNS
    -------
    resolution : string
      The coupling, currently N, NLS or ACX1 (returned in upper case)
    ihdu : int
      The HDU with the cross section data for the ion. Set to -1 for none.
    """

    ihdu = numpy.where( (self.crosssectiondata['INDEX'].data['Z']==Z) &\
                        (self.crosssectiondata['INDEX'].data['z1']==z1))[0]


    if len(ihdu) == 1:
      ihdu = ihdu[0]
      try:
        resolution = self.crosssectiondata['INDEX'].data['resn'][ihdu]
      except:
        resolution = self.crosssectiondata['INDEX'].data['resn'][ihdu].decode('ascii')

    else:
      ihdu = -1
      resolution = 'ACX1'

    if DEBUG:
      print("crosssectiondata type = %s, hdu = %i"%(resolution, ihdu))
    return resolution.upper(), ihdu


  def set_collisiontype(self, colltype, collunits='default'):
    """
    Set the collision type and units

    PARAMETERS
    ----------
    colltype : int
      Parameter for provided collision type
      Collision type - 1=energy/mass of center of mass
      Collision type - 2=velocity of center of mass
      Collision type - 3=velocity of donor
      Collision type - 4=velocity of receiver

    collunits : string, optional
      Units of collision paramter. Defaults to "kev/u" for colltype=1, "cm/s" for others


    RETURNS
    -------
    None
    """
    if colltype == 1:
      if collunits == 'default':
        self.collisiontype = 1
        self.collisionunits = 'kev/u'
      else:
        if collunits.lower() in ['kev/amu','kev/u','ev/u','ev/amu']:
          self.collisiontype = 1
          self.collisionunits = collunits.lower()
        else:
          print("Error: unknown units %s for center of mass collision energy")

    elif colltype in [2,3,4]:
      if collunits == 'default':
        self.collisiontype = colltype
        self.collisionunits = 'cm/s'
      else:
        if collunits.lower() in ['cm/s']:
          self.collisiontype = colltype
          self.collisionunits = collunits.lower()
        elif collunits.lower() in ['km/s']:
          self.collisiontype = colltype
          self.collisionunits = collunits.lower()
        else:
          print("Error: unknown units %s for center of mass collision energy")

    else:
      print("Error: unknown collision type: ", colltype,", should be 1, 2, 3 or 4")

  def set_collisionparam(self, collisionparam):
    """
    Set the collision velocity or energy

    PARAMETERS
    ----------
    collisionparam : float
      The collision velocity or energy. Units and parameter type are set in ACXModel.set_collisiontype

    RETURNS
    -------
    None
    """

    # to store the center of mass parameters
    self.collenergy={} # C.o.M. energies in kev/u
    self.collvelocity={} # C.o.M. velocities in cm/s

    for Z in self.elements:
      if self.collisiontype==1: # Energy of center of mass provided
        if self.collisionunits.lower() in ['kev/amu', 'kev/u']:
          self.collenergy[Z] = collisionparam*1.0
          self.collvelocity[Z] = 1e5* numpy.sqrt(4786031.3*self.collenergy[Z]/25.)
        elif self.collisionunits.lower() in ['ev/u', 'ev/amu']:
          self.collenergy[Z] = collisionparam/1000.0
          self.collvelocity[Z] = 1e5* numpy.sqrt(4786031.3*self.collenergy[Z]/25.)
        else:
          print("set_collisionparam: unknown collision units %s"%(self.collisionunits))
          return

      elif self.collisiontype==2:
        if self.collisionunits.lower() == 'cm/s':
          # this is the reduced mass velocity
          self.collvelocity[Z] = collisionparam*1.0
        elif self.collisionunits.lower() == 'km/s':
          self.collvelocity[Z] = collisionparam*1e5
        self.collenergy[Z] = (25/4786031.3) * (self.collvelocity[Z]/1e5)**2

      elif self.collisiontype == 3:
          # donor ion is moving
        if self.collisionunits.lower() == 'cm/s':
          self.collvelocity[Z] = self.donormass*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
        elif self.collisionunits.lower() == 'km/s':
          self.collvelocity[Z] = 1e5*self.donormass*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
        self.collenergy[Z] = 25 * (self.collvelocity[Z]/1e5)**2/4786031.3

      elif self.collisiontype == 4:
          # receiver ion is moving
        if self.collisionunits.lower() == 'cm/s':
          self.collvelocity[Z] = pyatomdb.atomic.Z_to_mass(Z)*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
        elif self.collisionunits.lower() == 'km/s':
          self.collvelocity[Z] = 1e5*pyatomdb.atomic.Z_to_mass(Z)*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
        self.collenergy[Z] = 25 * (self.collvelocity[Z]/1e5)**2/4786031.3

      else:
        print("*** ERROR: Unknown collision unit %s: should be kev/amu, km/s or  cm/s ***" %(self.collisionunits))
        return


  def set_ebins(self, ebins, ebins_checksum = False):
    """
    Set the energy bins for the spectrum being returned.

    PARAMETERS
    ----------
    ebins : array(float)
      Energy bin edges (keV)
    ebins_checksum : string, optional
      The hex digest of the md5 sum of ebins. Used to check for changes.
    """

    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()

    self.ebins = ebins

    if ebins_checksum != self.ebins_checksum:
#      print("setting ebins")
      self.ebins_checksum = ebins_checksum
      for Z in self.spectra.keys():
        for z1 in self.spectra[Z].keys():


          self.spectra[Z][z1].set_ebins(self.ebins, ebins_checksum=self.ebins_checksum)


  def calc_crosssection(self, collparam):
    """
    Calculate the total and partial cross sections by n, l, S.

    PARAMETERS
    ----------

    collparam : float
      Collision energy, or velocity, as determined by set_collisiontype, in kev/u,  cm/s or km/s

    RETURNS
    -------
    TBD

    """



  def calc_spectrum(self, collparam, Tbroaden, vbroaden, ebroaden):
    """
    Calculate the spectrum we want

    PARAMETERS
    ----------

    collparam : float
      Collision energy, or velocity, as determined by set_collisiontype, in kev/u,  cm/s or km/s

    RETURNS
    -------
    self.emiss : array(float)
      The emissivity in photons cm^3 bin-1 s-1
    """

    # get the velocity etc
    self.set_collisionparam(collparam)
    self.initialize_CXIonspectrum()

    # set up return array for spectrum
    if DEBUG:
      self.emiss_debug = {}
    self.emiss = numpy.zeros(len(self.ebins)-1, dtype=float)
#    print("IONFRAC:", self.ionfrac)
    #print("elements: ", self.elements)
    for iZ, Z in enumerate(self.elements):
      if self.abund[Z] > 0.0:
        for z1 in range(2, Z+2):# z1 here is the recombin*ing* ion charge +1
          if self.recombtype == SINGLE_RECOMBINATION:
            ionf = self.ionfrac[Z][z1-1] #ionfrac is indexed from 0, not 1
#            print("Z=",Z,", z1-1=", z1-1, ", ionf=", ionf)
          elif ((self.recombtype == FULL_RECOMBINATION) | (self.recombtype == FULL_NORMALIZED_RECOMBINATION)):
            ionf = sum(self.ionfrac[Z][z1-1:])
          else:
            raise ValueError("Invalid recombtype ", self.recombtype)
          if self.abund[Z]*ionf > 1e-10:
#            print("Z=%i, z1=%i"%(Z, z1))

            # set the energy bins for the spectrum
            self.spectra[Z][z1].set_ebins(self.ebins, ebins_checksum=self.ebins_checksum)
            self.spectra[Z][z1].set_recombtype(self.recombtype)

            if DEBUG:
              if not Z in self.emiss_debug.keys():
                self.emiss_debug[Z] = {}

              self.emiss_debug[Z][z1] = self.spectra[Z][z1].calc_spectrum(self.collenergy[Z], self.collvelocity[Z], Tbroaden, vbroaden, ebroaden)

              self.emiss +=  self.emiss_debug[Z][z1] * self.abund[Z] * ionf

            else:

              self.emiss += self.spectra[Z][z1].calc_spectrum(self.collenergy[Z], self.collvelocity[Z], Tbroaden, vbroaden, ebroaden) *\
                     self.abund[Z] * ionf


    return self.emiss

  def initialize_CXIonspectrum(self):
    """
    Called to ensure that all the ion spectra exist and are ready to use
    """

    for iZ, Z in enumerate(self.elements):
      if self.abund[Z] > 0.0:
        if not Z in self.spectra.keys():
          self.spectra[Z]={}
        for z1 in range(2, Z+2):# z1 here is the recombin*ing* ion charge +1
          if self.recombtype == SINGLE_RECOMBINATION:
            ionf = self.ionfrac[Z][z1-1] #ionfrac is indexed from 0, not 1

          elif ((self.recombtype == FULL_RECOMBINATION) |\
                (self.recombtype == FULL_NORMALIZED_RECOMBINATION)):
            ionf = sum(self.ionfrac[Z][z1-1:])
          else:
            raise ValueError("Invalid recombtype ", self.recombtype)
          if self.abund[Z]*ionf > 1e-10:


            if not z1 in self.spectra[Z].keys():
              # Initialize new CXIonSpectrum object for this ion
              resolution, ihdu = self.find_crosssection_type(Z,z1)

              if resolution=='ACX1':
                self.spectra[Z][z1] = CXIonSpectrum_ACX1(Z,z1,ihdu, \
                  self.linedata, self.contdata,\
                  acxmodel = self.acxmodel,\
                  donor = self.donor,\
                  receivermass = pyatomdb.atomic.Z_to_mass(Z),\
                  donormass = self.crosssectiondata['INDEX'].header['DONMASS'])


              elif resolution=='N':
                self.spectra[Z][z1] = CXIonSpectrum_N(Z,z1, self.crosssectiondata[ihdu+2].data, \
                  self.linedata, self.contdata,\
                  donor = self.crosssectiondata['INDEX'].header['DONMASS'],\
                  receivermass = self.crosssectiondata[ihdu+2].header['RECMASS'],\
                  donormass = self.crosssectiondata['INDEX'].header['DONMASS'])

              elif resolution=='NLS':
                self.spectra[Z][z1] = CXIonSpectrum_NLS(Z,z1, self.crosssectiondata[ihdu+2].data, \
                  self.linedata, self.contdata,\
                  donor = self.crosssectiondata['INDEX'].header['DONMASS'],\
                  receivermass = self.crosssectiondata[ihdu+2].header['RECMASS'],\
                  donormass = self.crosssectiondata['INDEX'].header['DONMASS'])

  def calc_linelist(self, collparam, specrange, specunit='A', redshift=0.0):
    """
    Calculate the spectrum we want

    PARAMETERS
    ----------

    collparam : float
      Collision energy, or velocity, as determined by set_collisiontype, in kev/u,  cm/s or km/s
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)

    RETURNS
    -------
    linelist : numpy.array
      Numpy array of lines and emissivity, custom dtype
    """

    # get the velocity etc
    self.set_collisionparam(collparam)


    # set up return array for spectrum
#    if DEBUG:
#      self.emiss_debug = {}
#    self.emiss = numpy.zeros(len(self.ebins)-1, dtype=float)

    self.initialize_CXIonspectrum()

    # zero out the line list
    self.linelist = numpy.zeros(0, dtype = llistdtype)

    for iZ, Z in enumerate(self.elements):
      if self.abund[Z] > 0.0:
        for z1 in range(2, Z+2):# z1 here is the recombin*ing* ion charge +1
          if self.recombtype == SINGLE_RECOMBINATION:
            ionf = self.ionfrac[Z][z1-1] #ionfrac is indexed from 0, not 1

          elif ((self.recombtype == FULL_RECOMBINATION) |\
                (self.recombtype == FULL_NORMALIZED_RECOMBINATION)):
            ionf = sum(self.ionfrac[Z][z1-1:])
          else:
            raise ValueError("Invalid recombtype ", self.recombtype)
          if self.abund[Z]*ionf > 1e-10:

            self.spectra[Z][z1].set_recombtype(self.recombtype)

            tmp = self.spectra[Z][z1].calc_linelist(\
                                     self.collenergy[Z], \
                                     self.collvelocity[Z], \
                                     specrange, specunit=specunit, \
                                     redshift=redshift)
            tmp['EPSILON'] *= self.abund[Z] * ionf
            self.linelist=numpy.append(self.linelist,tmp)



    return self.linelist


  def calc_line_emissivity(self, Z, z1_in, up, lo, collvalue=None):
    """
    Calculate the emissivity of a specific line.

    PARAMETERS
    ----------
    collvalue : float
      The collision parameter (kev/u or cm/s, depending) to calculate the spectrum
    Z : int
      element charge
    z1 : int
      recombining ion charge +1.
    up : int
      upper level of transition
    lo : int
      lower level of transition

    RETURNS
    -------
    ret : array(float)
      The spectrum, in ph cm^3 s-1 bin-1
    """

    if collvalue is not None:
      self.set_collisionparam(collvalue)

    emissivity = 0.0 # default number
    z1_ion = z1_in +1 # If you want a specific line, say He-like 7 to 1, you need
                   # to start with the ion  abundance of the H-like so it can recombine
    ret = {'Lambda': 0.0,
           'Epsilon': 0.0,
           'up': up,
           'lo': lo}

    self.initialize_CXIonspectrum()

    if self.abund[Z] > 0.0:
      if not Z in self.spectra.keys():
        self.spectra[Z]={}
      for z1 in [z1_ion]:# z1 here is the recombin*ing* ion charge +1
        if self.recombtype == SINGLE_RECOMBINATION:
          ionf = self.ionfrac[Z][z1-1]
        elif self.recombtype == FULL_RECOMBINATION:
          ionf = sum(self.ionfrac[Z][z1-1:])
        else:
          raise ValueError("Invalid recombtype ", self.recombtype)
        if self.abund[Z]*ionf > 1e-40: #ionfrac is indexed from 0, not 1
#          print("Z=%i, z1=%i"%(Z, z1))
            # set the energy bins for the spectrum
          self.spectra[Z][z1].set_ebins(self.ebins, ebins_checksum=self.ebins_checksum)

          if DEBUG:
            if not Z in self.emiss_debug.keys():
              self.emiss_debug[Z] = {}

#              self.emiss_debug[Z][z1] = self.spectra[Z][z1].calc_spectrum(self.ebins, self.collenergy[Z], self.collvelocity[Z], self.linedata, self.contdata, self.acxmodel, Tbroaden, vbroaden)
            self.emiss_debug[Z][z1] = self.spectra[Z][z1].calc_spectrum(self.collenergy[Z], self.collvelocity[Z], Tbroaden, vbroaden, ebroaden)

            self.emiss +=  self.emiss_debug[Z][z1] * self.abund[Z] * ionf

          else:
            line_emissivity = self.spectra[Z][z1].calc_line_emissivity(self.collenergy[Z], self.collvelocity[Z], up, lo)
            ret['Epsilon'] +=line_emissivity['Epsilon']* self.abund[Z] * ionf
            if line_emissivity['Lambda']>0:
              ret['Lambda']=line_emissivity['Lambda']

    return(ret)





  def set_temperature(self, temperature):
    """
    Recalculate the ionization balance based on equilibrium electron temperature

    PARAMETERS
    ----------
    temperature: float
      electron temperature in keV

    RETURNS
    -------
    None

    """

    # set the numbersset_collisionparam
    self.ionfrac_from_temperature = True
    self.temperature = temperature

    # calculate ionization balance
    self.calc_ionfrac_equilibrium()

  def set_ionfrac(self, ionfrac):
    """
    Recalculate the ionization balance based on equilibrium electron temperature

    PARAMETERS
    ----------
    ionfrac: dict of arrays
      ionization fraction, e.g. ionfrac[8]=numpy.array([0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.4, 0.2, 0.0])

    RETURNS
    -------
    None

    """

    # set the numbersset_collisionparam
    self.ionfrac_from_temperature = False
    self.temperature = False

    # calculate ionization balance
    self.ionfrac = ionfrac


  def set_acxmodel(self, acxmodel):
    """
    Set the ACX spectrum type

    PARAMETERS
    ----------
    acxmodel : int
      The acxmodel (between 1 and 8)

    """
    self.acxmodel=acxmodel

    for Z in self.spectra.keys():
      for z1 in self.spectra[Z].keys():
        self.spectra[Z][z1].set_acxmodel(self.acxmodel)


  def set_recombtype(self, recombtype):
    """
    Set the recombination type

    PARAMETERS
    ----------
    recombtype : int
      The type of recombination (1=single, 2=all the way to neutral,
      3 = same as 2, but renomalized by total cross section)

    """
    self.recombtype = recombtype

  def set_abund(self, abund, elements=None):
    """
    Set the elemental abundance, also in each donor model

    PARAMETERS
    ----------
    abund : array(float)
      Abundances, relative to defulat

    RETURNS
    -------
    None
    """

    try:
      if elements==None:
        elements = self.elements
      else:
        elements, elissvec = pyatomdb.util.make_vec(elements)
    except ValueError:
      elements, elissvec = pyatomdb.util.make_vec(elements)

    abundvec, aisvec = pyatomdb.util.make_vec(abund)

    if len(abundvec) != len(elements):
      if len(abundvec)==1:
        abundvec = abundvec[0]*numpy.ones(len(elements))

      else:
        raise(ValueError,"Error: specified incompatible number of elements from set_abund")
        return

    for i in range(len(abundvec)):
      self.abund[elements[i]] = abundvec[i]



class CXIonSpectrum():
  """
  Class to store and prepare each ion's spectrum. Will have a
  subset of CXShellSpectrum for each shell.
  Can provide it a set of energy bins and get a spectrum back.
  This will store all the emissivity data.

  """
  import pyatomdb, numpy, os, hashlib

  def __init__(self, Z, z1, ebins, crosssectiondata,\
               linedata, contdata, cxdata):
    """
    initiatlize

    PARAMETERS
    ----------0.1
    Z : int
      Nuclear Charge
    z1 : int
      Charge +1 of *recombining* ion
#    collision_energy : float
#      Collision energy in keV/u
#    ebins : array(float)
#      Energy bins.
    crossssectiondata : HDUList
      Cross section data for all ions.
    """

    # set some values
    self.Z = Z
    self.z1 = z1
    self.donor = donor
    self.Ip_donor = linedata[1].header['Ip_D']

    # find the matching cross section data
    ihdu = numpy.where( (crosssectiondata['INDEX'].data['Z']==self.Z) &\
                        (crosssectiondata['INDEX'].data['z1']==self.z1))[0]
    if len(ihdu) == 1:
      self.crosssectiondata = crosssectiondata[ihdu+2].data
      self.coupling = crosssectiondata['INDEX'].data['resn'][i].decode('ascii')
      self.DonorMass = crosssectiondata[ihdu+2].header['DONMASS']
      self.RecvMass = crosssectiondata[ihdu+2].header['RECMASS']

    else:
      self.coupling='ACX1'
      self.RecvMass=pyatomdb.atomic.Z_to_mass(self.Z)
      self.DonorMass=crosssectiondata['INDEX'].header['DONMASS']

    # now find the matching parts of emisssivity
    self.linedataindexes = {}
    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0]

    for ii in i:
      n = linedata[1].data[ii]['n']
      if not n in self.linedataindexes.keys():
        self.linedataindexes[n] = {}
      l = linedata[1].data[ii]['l']


      if self.coupling=='nls':
        if not l in self.linedataindexes[n].keys():
          self.linedataindexes[n][l]={}
        s = linedata[1].data[ii]['s2p1']
        if not s in self.linedataindexes[n][l].keys():
          self.linedataindexes[n][l][s]=ii+2

      else:
        if not l in self.linedataindexes[n].keys():
          self.linedataindexes[n][l]=ii+2

    # now find the matching parts of emisssivity
    self.contdataindexes = {}
    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0]

    for ii in i:
      n = contdata[1].data[ii]['n']
      if not n in self.contdataindexes.keys():
        self.contdataindexes[n] = {}
      l = contdata[1].data[ii]['l']


      if self.coupling=='nls':
        if not l in self.contdataindexes[n].keys():
          self.contdataindexes[n][l]={}
        s = contdata[1].data[ii]['s2p1']
        if not s in self.contdataindexes[n][l].keys():
          self.contdataindexes[n][l][s]=ii+2

      else:
        if not l in self.contdataindexes[n].keys():
          self.contdataindexes[n][l]=ii+2




      if not contdata[1].data[ii]['n'] in self.contdataindexes.keys():
        self.contdataindexes[contdata[1].data[ii]['n']] = {}



  def set_ebins(self, ebins, ebins_checksum=False):
    """
    Set the energy bins for the spectrum being returned.

    PARAMETERS
    ----------
    ebins : array(float)
      Energy bin edges (keV)
    ebins_checksum : string, optional
      The hex digest of the md5 sum of ebins. Used to check for changes.
    """

    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()


    if not ebins_checksum == self.ebins_checksum:
      # need to reset everything

      self.ebins = ebins # tracks for changes to ebins
      self.ebins_checksum = ebins_checksum # tracks for changes to ebins

      self.reset_ebins(ebins, ebins_checksum=self.ebins_checksum)


  def return_spectrum(self, ebins, collision_energy, linedata, contdata,\
                      ebins_checksum=False):
    """
    Return the spectrum on the energy bins

    """

    # For each and every data set, calculate the relevant cross sections

    if self.couping=='ACX1':
      # We have old style ACX coupling
      q = 1.0*(self.z1-1)
      Ip_d =self.Ip_donor
      nprime = q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                        (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))

      l = -1 * (self.acxmodel%4)

      # calculate n shell distribution
      if 1 <= self.acxmodel <=4:

        # need to do just 1 n shell, the closest.
        n = [int(numpy.round(nprime))]
        nfrac=[1.0]

      else:
        n = [int(numpy.floor(nprime)), int(numpy.ceil(nprime))]
        nfrac = [1-(nprime%1), nprime%1]


      for inn, nn in enumerate(n):
        # check if this already exists.

        ####HERE####
        calc_single_spectrum(ebins, nn, l,linedata, contdata)

  def set_acxmodel(self, acxmodel):
    """
    Set the ACX spectrum type

    PARAMETERS
    ----------
    acxmodel : int
      The acxmodel (between 1 and 8)

    """
    self.acxmodel=acxmodel


  #copy expand_E_grid into CXIonSpectrum
  def expand_E_grid(self, eedges, Econt_in_full, cont_in_full):

    """
    Code to expand the compressed continuum onto a series of bins.

    Parameters
    ----------
    eedges : float(array)
      The bin edges for the spectrum to be calculated on, in units of keV
    Econt_in_full: float(array)
      The compressed continuum energies (keV)
    cont_in_full: float(array)
      The compressed continuum emissivities (ph cm3 s-1 keV-1)

    Returns
    -------
    float(array)
      len(bins)-1 array of continuum emission, in units of \
      photons cm^3 s^-1 bin^-1
    """

  #  History
  #  -------
  #  Version 0.1 - initial release
  #    Adam Foster July 17th 2015

    import scipy.integrate
    cont_in = cont_in_full
    Econt_in = Econt_in_full
    n = len(cont_in_full)

    # Append desired bin edges to the end of this array
    E_all = numpy.append(Econt_in, eedges)

    # interpolate the continuum at these bin edges
    cont_tmp = numpy.interp(eedges, Econt_in, cont_in)

    # add continuum at bin edges to the end
    C_all = numpy.append(cont_in, cont_tmp)

    # get sorted indexes for energy
    iord = numpy.argsort(E_all)

    # reorder the arrays in energy order
    E_all = E_all[iord]
    C_all = C_all[iord]

    # pull out the points which are from the output grid (iord > npoints)
    ihi = numpy.where(iord>=n)[0]

    # do cumulative sum integration
    cum_cont = scipy.integrate.cumulative_trapezoid(C_all, E_all, initial=0)

    # extract the cumulative sum at the output points
    C_out = cum_cont[ihi]

    # convert to flux per bin
    cont = C_out[1:]-C_out[:-1]

    return cont


  def set_recombtype(self, recombtype):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """

    if not recombtype in [SINGLE_RECOMBINATION,FULL_RECOMBINATION, FULL_NORMALIZED_RECOMBINATION]:
      print("ERROR: invalid recombtype ", recombtype)

    else:
      self.recombtype = recombtype



class CXIonSpectrum_ACX1(CXIonSpectrum):
  """
  This is a class for ACX1 model data

  ATTRIBUTES
  ----------
  ebins : array(float)
    The energy bin in keV
  ebins_m5sum : md5hash
    hex digest of ebins

  """
  import pyatomdb, numpy, os, hashlib

  def __init__(self, Z,z1,crosssectionhdu, \
               linedata,\
               contdata,\
               acxmodel = 8,\
               donor = False,\
               receivermass = False,\
               donormass = False):
    """
    initialize with some important variaset_collisionparambles

    Z : int
      recombining ion nuclear charge
    z1 : int
      recombining ion charge +1
    crosssectionhdu: int
      The number (indexed from 0) of the HDU with the cross section data
    linedata : hdulist
      The line emissivity data for this ion
    contdata : hdulist
      The continuum emissivity data for this ion
    acxmodel : int, optional (default = 8)
      The ACX model number in use.
    donor : string
      The donor symbol (e.g. H, He, or H2O)
    receivermass : float
      Mass in a.m.u. of the recombining particle
    donormass : float
      Mass in a.m.u. of the donor particle
    """

    # store all the info

    self.Z = Z
    self.z1 = z1

    self.Ip_donor = linedata[1].header['Ip_D']
    self.crosssectionhdu = crosssectionhdu
    self.donor = donor
    self.receivermass = receivermass
    self.donormass = donormass
    self.acxmodel = acxmodel

    # set up some empty variables, to be fleshed out
    self.ebins_checksum = "0" # store hash of ebins
    # find, and store, the relevant HDUs from the line and continuum datafiles.


    #find the hdu location for the matching Z and z1 for line and continuum
    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0]

    if len(i) > 0:
      self.ionlinedata = linedata[i[0]+2].data
    else:
      self.ionlinedata = numpy.zeros(0, dtype = linedata[3].data.dtype)



    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0]

    if len(i) > 0:
      self.ioncontdata = contdata[i[0]+2].data

    #create arrays to hold the n l s values
      self.n = numpy.zeros(len(self.ioncontdata))
      self.l = numpy.zeros(len(self.ioncontdata))
      self.s = numpy.zeros(len(self.ioncontdata))

    #find the n l s values from the header and write to arrays
      for j,item in enumerate(self.ioncontdata):
         headstr = "CONT_*"
         headstr = headstr.replace("*", str(j))
         nlsstr = contdata[i[0]+2].header[headstr]
         nlslist = nlsstr.replace('n',' ').replace('l',' ').replace('s',' ').split()

         self.n[j] = nlslist[0]
         self.l[j] = nlslist[1]
         self.s[j] = nlslist[2]

    else:
      self.ioncontdata = cocodata = numpy.zeros(0,dtype=pyatomdb.apec.generate_datatypes('continuum', ncontinuum=0, npseudo=0))
      self.n = numpy.zeros(0)
      self.l = numpy.zeros(0)
      self.s = numpy.zeros(0)

    #to store the location of all the HDU data
    self.linedata = linedata
    self.contdata = contdata

    self.spectra ={}

    self.ebins=0.0
    self.ebins_checksum="0"

  def set_ebins(self, ebins, ebins_checksum=False):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """



    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()
    if ebins_checksum != self.ebins_checksum:
      self.ebins = ebins
      self.ebins_checksum = ebins_checksum

    for n in self.spectra.keys():
      for l in self.spectra[n].keys():

        self.spectra[n][l].set_ebins(ebins, ebins_checksum= self.ebins_checksum)



#  def set_acxmodel(self, acxmodel):
#    """
#    Set the ACX spectrum typ
#
#    PARAMETERS
#    ----------
#    acxmodel : int
#      The acxmodel (between 1 and 8)
#
#    """
#    self.set_acxmodel=acxmodel

#  def calc_spectrum(self, ebins, collenergy, collvelocity, linedata, contdata, acxmodel, Tbroaden, vbroaden):
  def calc_crosssection(self, collenergy):
    """
    Calculate cross sections
    """
#    +=====+===================+==================================+
#    |Value| n distribution    | l, L distribution                |
#    +-----+-------------------+----------------------------------+                                                                               |
#    |  1  | one n shell       | even distribution by  l.         |
#    |  2  | one n shell       | statistical distribution by l.   |
#    |  3  | one n shell       | Landau-Zener distribution by  l. |
#    |  4  | one n shell       | Separable distribution by l.     |
#    |  5  | weighted 2 shells | even distribution by  l.         |
#    |  6  | weighted 2 shells | statistical distribution by l.   |
#    |  7  | weighted 2 shells | Landau-Zener distribution by l.  |
#    |  8  | weighted 2 shells | Separable distribution by l.     |
#    +-----+-------------------+----------------------------------+
    q = 1.0*(self.z1-1)
    Ip_d =self.Ip_donor
    nprime = q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                  (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))


    l = -1 * (self.acxmodel%4)

    # if l = -4, this returns 0. Correct back to -4
    if l==0: l=-4

    if 1 <= self.acxmodel <=4:

        # need to do just 1 n shell, the closest.
      n = [int(numpy.round(nprime))]
      nfrac=[1.0]
    else:
      n = [int(numpy.floor(nprime)), int(numpy.ceil(nprime))]
      nfrac = [1-(nprime%1), nprime%1]

    #create array to hold nfrac in same slot as choosen n l s
    capture_frac = numpy.zeros(len(self.n), dtype = float)

    for b in range(len(n)):
      for k,item in enumerate(self.n):
        if self.n[k] == n[b] and self.l[k] == l:
          capture_frac[k] = nfrac[b]



    ret = {}
    ret['n'] = self.n
    ret['l'] = self.l
    ret['S'] = numpy.zeros(len(self.n), dtype=int)
    ret['S'][:] = -1
    ret['sigma'] = capture_frac*UNIVERSAL_CX_CROSSSECTION
    ret['sigma_total'] = sum(capture_frac) # should be UNIVERSAL_CX_CROSSSECTION

    return(ret)

  def calc_spectrum(self, collenergy, collvelocity, Tbroaden, vbroaden, ebroaden):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    linedata : hdulist
      The line emissivity data
    contdata : hdulist
      The continuum emissivity data
    acxmodel : int
      n, l shell distribution, between 1 and 8:


    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """
#    +=====+===================+==================================+
#    |Value| n distribution    | l, L distribution                |
#    +-----+-------------------+----------------------------------+                                                                               |
#    |  1  | one n shell       | even distribution by  l.         |
#    |  2  | one n shell       | statistical distribution by l.   |
#    |  3  | one n shell       | Landau-Zener distribution by  l. |
#    |  4  | one n shell       | Separable distribution by l.     |
#    |  5  | weighted 2 shells | even distribution by  l.         |
#    |  6  | weighted 2 shells | statistical distribution by l.   |
#    |  7  | weighted 2 shells | Landau-Zener distribution by l.  |
#    |  8  | weighted 2 shells | Separable distribution by l.     |
#    +-----+-------------------+----------------------------------+

    # This will be returned NOT accounting for ionization balance or
    # elemental abundance. It is assumed that this will be handled at
    # a higher level, and indeed, since this call is quite time consuming,
    # it shouldn't be made at all for unnecessary ions.


    # job is to create a single unified return.

    # q = 1.0*(self.z1-1)
    # Ip_d =self.Ip_donor
    # nprime = q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                  # (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))

    # l = -1 * (self.acxmodel%4)
    # # if l = -4, this returns 0. Correct
    # if l==0: l=-4

      # # calculate n shell distribution
    # if 1 <= self.acxmodel <=4:

        # # need to do just 1 n shell, the closest.
      # n = [int(numpy.round(nprime))]
      # nfrac=[1.0]
    # else:
      # n = [int(numpy.floor(nprime)), int(numpy.ceil(nprime))]
      # nfrac = [1-(nprime%1), nprime%1]

    # #create array to hold nfrac in same slot as choosen n l s
    # capture_frac = numpy.zeros(len(self.n), dtype = float)
# #    print('n',n, 'nfrac:', nfrac)
    # for b in range(len(n)):
      # for k,item in enumerate(self.n):
        # if self.n[k] == n[b] and self.l[k] == l:
          # capture_frac[k] = nfrac[b]

# #    print('capture_frac:', capture_frac)


    cross_section = self.calc_crosssection(collenergy)
    Coutarraynl = cross_section['sigma']

    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    #create array to return spectrum

    #######line broadening
    self.broaden_limit = 1e-30
    self.broaden_object = _Gaussian_CDF()
    spec = numpy.zeros(len(self.ebins)-1, dtype=float)

    if ((Tbroaden <= 0.0) & \
          (vbroaden <=0.0) &\
          (ebroaden <=0.0)):

      if len(self.ionlinedata) > 0:
        #turn line wavelengths into energies
        en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
        #puts lines into the right spectral bins using epsilon array
        a,b= numpy.histogram(en, bins = self.ebins, weights = new_epsilon) #new_epsilon
        #add the line emission per bin to the spectrum
        spec += a

    elif len(self.ionlinedata) > 0:

      en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
      ind = new_epsilon>self.broaden_limit
      nonind = ~ind

      # calculate the widths of the strong lines
      llist = self.ionlinedata[ind]

      # get a raw dictionary of masses in amu
      masslist = pyatomdb.atomic.Z_to_mass(1,raw=True)

      #assign values of temperature and velocity broadening
      T = Tbroaden
      Tb = pyatomdb.util.convert_temp(T, 'K','keV')*pyatomdb.const.ERG_KEV/(masslist[llist['Element']]*1e3*pyatomdb.const.AMUKG)

      if vbroaden <0:
        vbroaden = 0.0
        vb=0.0
      else:
        vb = (vbroaden * 1e5)**2

      if ebroaden <0:
        ebroaden = 0.0
        eb=0.0
      else:
        eb = ebroaden**2

      wcoeff = numpy.sqrt(Tb+vb+eb) / (pyatomdb.const.LIGHTSPEED*1e2)

      elines = en[ind]
      n_e = new_epsilon[ind]
      width = wcoeff*elines

      # Filter out lines more than NSIGMALIMIT sigma outside the range
      NSIGMALIMIT=4
      eplu = elines+NSIGMALIMIT*width
      eneg = elines-NSIGMALIMIT*width
      emax = max(self.ebins)
      emin = min(self.ebins)
      # identify all the good lines!
      igood = numpy.where(((elines >= emin) & (eneg < emax))  |\
                ((elines < emax) & (eplu > emin)))[0]

      tmpspec = numpy.zeros(len(self.ebins))
#      t0 = time.time()
      for iline in igood:
        #if width
#        tmpspec += self.broaden_object.broaden(pyatomdb.const.HC_IN_KEV_A/llist['Lambda'][iline],\
#                         width[iline],self.ebins)*new_epsilon[iline]
        tmpspec += self.broaden_object.broaden(elines[iline],\
                         width[iline],self.ebins)*n_e[iline]
#      t1 = time.time()
#      print("  Time spent broadening: %.3f"%(t1-t0))
      spec=tmpspec[1:]-tmpspec[:-1]
      s,z = numpy.histogram(en[nonind], \
                            bins = self.ebins,\
                            weights = new_epsilon[nonind])
      spec+=s

    #add to spectrum for the capture_frac n l s values
    if len(self.ioncontdata) > 0:
      for p,item in enumerate(self.ioncontdata):
        if Coutarraynl[p] > 0.0:

          if self.ioncontdata['N_cont'][p] > 2:
            ncont = self.ioncontdata['N_Cont'][p] #shorten array
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Cont'][p][:ncont], self.ioncontdata['Continuum'][p][:ncont]) * (Coutarraynl[p])
          if self.ioncontdata['N_Pseudo'][p] > 2:
            npse = self.ioncontdata['N_Pseudo'][p]
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Pseudo'][p][:npse], self.ioncontdata['Pseudo'][p][:npse]) * (Coutarraynl[p])

# This line for renormalization
#    if self.recombtype != FULL_NORMALIZED_RECOMBINATION:
#      spec = spec * UNIVERSAL_CX_CROSSSECTION
#    else:
#      print("rescaling spectrum by %e"%(UNIVERSAL_CX_CROSSSECTION))

    self.ebins_checksum = self.ebins_checksum
    self.spectrum = spec
    self.spectrum_ready = True

    return self.spectrum * collvelocity

  def calc_linelist(self, collenergy, collvelocity, specrange, specunit='A', redshift=0.0):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)




    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """
#    +=====+===================+==================================+
#    |Value| n distribution    | l, L distribution                |
#    +-----+-------------------+----------------------------------+                                                                               |
#    |  1  | one n shell       | even distribution by  l.         |
#    |  2  | one n shell       | statistical distribution by l.   |
#    |  3  | one n shell       | Landau-Zener distribution by  l. |
#    |  4  | one n shell       | Separable distribution by l.     |
#    |  5  | weighted 2 shells | even distribution by  l.         |
#    |  6  | weighted 2 shells | statistical distribution by l.   |
#    |  7  | weighted 2 shells | Landau-Zener distribution by l.  |
#    |  8  | weighted 2 shells | Separable distribution by l.     |
#    +-----+-------------------+----------------------------------+

    # convert input wavelength or energy range to angstroms
    wave = pyatomdb.util.convert_spec(specrange, specunit, 'A')

    cross_section = self.calc_crosssection(collenergy)
    Coutarraynl = cross_section['sigma']

    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    #create array to return spectrum
    l = numpy.zeros(len(self.ionlinedata), dtype=llistdtype)

    if len(self.ionlinedata) > 0:

      l['EPSILON'] = new_epsilon * collvelocity
      l['LAMBDA'] =  self.ionlinedata['LAMBDA']
      l['ENERGY'] = pyatomdb.const.HC_IN_KEV_A/l['LAMBDA']
      l['ELEMENT'] = self.ionlinedata['ELEMENT']
      l['ION'] = self.ionlinedata['ION']
      l['ION_DRV'] = self.ionlinedata['ION_DRV']
      l['UPPERLEV'] = self.ionlinedata['UPPERLEV']
      l['LOWERLEV'] = self.ionlinedata['LOWERLEV']

    self.linelist = l[(l['LAMBDA']>=wave[0]) & (l['LAMBDA']<=wave[1])]
    return self.linelist



  def calc_line_emissivity(self, collenergy, collvelocity, up, lo):

    emissivity = 0.0
    wv_out = 0.0
#    +=====+===================+==================================+
#    |Value| n distribution    | l, L distribution                |
#    +-----+-------------------+----------------------------------+                                                                               |
#    |  1  | one n shell       | even distribution by  l.         |
#    |  2  | one n shell       | statistical distribution by l.   |
#    |  3  | one n shell       | Landau-Zener distribution by  l. |
#    |  4  | one n shell       | Separable distribution by l.     |
#    |  5  | weighted 2 shells | even distribution by  l.         |
#    |  6  | weighted 2 shells | statistical distribution by l.   |
#    |  7  | weighted 2 shells | Landau-Zener distribution by l.  |
#    |  8  | weighted 2 shells | Separable distribution by l.     |
#    +-----+-------------------+----------------------------------+

    # This will be returned NOT accounting for ionization balance or
    # elemental abundance. It is assumed that this will be handled at
    # a higher level, and indeed, since this call is quite time consuming,
    # it shouldn't be made at all for unnecessary ions.


    # job is to create a single unified return.

    q = 1.0*(self.z1-1)
    Ip_d =self.Ip_donor
    nprime = q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                  (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))

    l = -1 * (self.acxmodel%4)
    # if l = -4, this returns 0. Correct
    if l==0: l=-4

      # calculate n shell distribution
    if 1 <= self.acxmodel <=4:

        # need to do just 1 n shell, the closest.
      n = [int(numpy.round(nprime))]
      nfrac=[1.0]
    else:
      n = [int(numpy.floor(nprime)), int(numpy.ceil(nprime))]
      nfrac = [1-(nprime%1), nprime%1]

    #create array to hold nfrac in same slot as choosen n l s
    capture_frac = numpy.zeros(len(self.n), dtype = float)

    for b in range(len(n)):
      for k,item in enumerate(self.n):
        if self.n[k] == n[b] and self.l[k] == l:
          capture_frac[k] = nfrac[b]

    if len(capture_frac)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*capture_frac,1)

    i = numpy.where((self.ionlinedata['UpperLev'] == up) &
                    (self.ionlinedata['LowerLev'] == lo))[0]
    if len(i)==0:
      out_wv = 0.0
      emissivity = 0.0
    else:
      out_wv = self.ionlinedata['Lambda'][i[0]]
      emissivity = new_epsilon[i[0]]* UNIVERSAL_CX_CROSSSECTION

    emissivity *= collvelocity

    ret = {'Lambda': out_wv,
           'Epsilon': emissivity,
           'up': up,
           'lo': lo}

    return ret



class CXIonSpectrum_NLS(CXIonSpectrum):
  """
  This is a class for nls resolved Kronos data
  """
  import pyatomdb, numpy, os, hashlib


  def __init__(self, Z,z1,crosssectiondata, \
               linedata,\
               contdata,\
               donor = False,\
               receivermass = False,\
               donormass = False):
    """
    initialize with some important variaset_collisionparambles

    Z : int
      recombining ion nuclear charge
    z1 : int
      recombining ion charge +1
    crosssectiondata: hdulist
      The  HDU with the cross section data
    linedata : hdulist
      The line emissivity data for this ion
    contdata : hdulist
      The continuum emissivity data for this ion
    acxmodel : int, optional (default = 8)
      The ACX model number in use.
    donor : string
      The donor symbol (e.g. H, He, or H2O)
    receivermass : float
      Mass in a.m.u. of the recombining particle
    donormass : float
      Mass in a.m.u. of the donor particle
    """

    # store all the info

    self.Z = Z
    self.z1 = z1

    self.Ip_donor = linedata[1].header['Ip_D']
    self.crosssectiondata = crosssectiondata
    self.donor = donor
    self.receivermass = receivermass
    self.donormass = donormass

    # set up some empty variables, to be fleshed out
    self.ebins_checksum = "0" # store hash of ebins
    # find, and store, the relevant HDUs from the line and continuum datafiles.

    #find the hdu location for the matching Z and z1 for line and continuum
    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0][0]

    self.ionlinedata = linedata[i+2].data

    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0][0]

    self.ioncontdata = contdata[i+2].data

    #find the n l s values from the header and write to arrays
    self.n = numpy.zeros(len(self.ioncontdata))
    self.l = numpy.zeros(len(self.ioncontdata))
    self.s = numpy.zeros(len(self.ioncontdata))

    for j,item in enumerate(self.ioncontdata):

       headstr = "CONT_*"
       headstr = headstr.replace("*", str(j))
       nlsstr = contdata[i+2].header[headstr]
       nlslist = nlsstr.replace('n',' ').replace('l',' ').replace('s',' ').split()

       self.n[j] = nlslist[0]
       self.l[j] = nlslist[1]
       self.s[j] = nlslist[2]

    #to store the location of all the HDU data
    self.linedata=linedata
    self.contdata=contdata


    self.spectra ={}
    self.ebins=0.0
    self.ebins_checksum="0"



  def calc_spectrum(self, collenergy, collvelocity, Tbroaden, vbroaden, ebroaden):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    linedata : hdulist
      The line emissivity data
    contdata : hdulist
      The continuum emissivity data
    acxmodel : int
      n, l shell distribution, between 1 and 8


    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """


    # This will be returned NOT accounting for ionization balance or
    # elemental abundance. It is assumed that this will be handled at
    # a higher level, and indeed, since this call is quite time consuming,
    # it shouldn't be made at all for unnecessary ions.


    # cycle through every single entry in cross section file
    #sigmadata = self.crosssectiondata

    #find energy range
    #emin = min(sigmadata['E_min'])
    #emax = min(sigmadata['E_max'])
#    print("NLS Spectrum calc")

    #create array that goes through the sigma files for the known distribution of n l s
    cross_section = self.calc_crosssection(collenergy)

    Coutarraynl = cross_section['sigma']
    #for every transition, multiplies the capture_frac with the emission, and sums the value. total emissivity
#    new_epsilon = numpy.sum(Coutarraynl * self.ionlinedata['Epsilon'], 1)
    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    #set up array to return spectrum
    spec = numpy.zeros(len(self.ebins)-1, dtype=float)

    #######line broadening
    self.broaden_limit = 1e-20
    self.broaden_object = _Gaussian_CDF()

    if ((Tbroaden <= 0.0) & \
          (vbroaden <=0.0)):

      if len(self.ionlinedata) > 0:
        #turn line wavelengths into energies
        en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
        #puts lines into the right spectral bins using epsilon array
        a,b= numpy.histogram(en, bins = self.ebins, weights = new_epsilon)
        #add the line emission per bin to the spectrum
        spec += a

    elif len(self.ionlinedata) > 0:

      en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
      ind = new_epsilon>self.broaden_limit
      nonind = ~ind

      # calculate the widths of the strong lines
      llist = self.ionlinedata[ind]

      # get a raw dictionary of masses in amu
      masslist = pyatomdb.atomic.Z_to_mass(1,raw=True)

      #assign values of temperature and velocity broadening
      T = Tbroaden
      Tb = pyatomdb.util.convert_temp(T, 'K','keV')*pyatomdb.const.ERG_KEV/(masslist[llist['Element']]*1e3*pyatomdb.const.AMUKG)

      if vbroaden <0:
        vbroaden = 0.0
        vb=0.0
      else:
        vb = (vbroaden * 1e5)**2

      wcoeff = numpy.sqrt(Tb+vb) / (pyatomdb.const.LIGHTSPEED*1e2)

      elines = en[ind]
      n_e = new_epsilon[ind]
      width = wcoeff*elines

      # Filter out lines more than NSIGMALIMIT sigma outside the range
      NSIGMALIMIT=4
      eplu = elines+NSIGMALIMIT*width
      eneg = elines-NSIGMALIMIT*width
      emax = max(self.ebins)
      emin = min(self.ebins)
      # identify all the good lines!

      igood = numpy.where(((elines >= emin) & (eneg < emax))  |\
                ((elines < emax) & (eplu > emin)))[0]

      spec = numpy.zeros(len(self.ebins))

      for iline in igood:
        #if width
#        spec += self.broaden_object.broaden(pyatomdb.const.HC_IN_KEV_A/llist['Lambda'][iline],\
#                         width[iline],self.ebins)*new_epsilon[iline]
        spec += self.broaden_object.broaden(elines[iline],\
                         width[iline],self.ebins)*n_e[iline]

      spec=spec[1:]-spec[:-1]

      s,z = numpy.histogram(en[nonind], \
                            bins = self.ebins,\
                            weights = new_epsilon[nonind])
      spec+=s
    #######

    #add to spectrum the continuum data for the Coutarraynl n l s values
    if len(self.ioncontdata) > 0:
      for p,item in enumerate(self.ioncontdata):
        if Coutarraynl[p] > 0.0:
          if self.ioncontdata['N_cont'][p] > 2:
            ncont = self.ioncontdata['N_Cont'][p] #shorten array
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Cont'][p][:ncont], self.ioncontdata['Continuum'][p][:ncont]) * (Coutarraynl[p])

          if self.ioncontdata['N_Pseudo'][p] > 2:
            npse = self.ioncontdata['N_Pseudo'][p]
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Pseudo'][p][:npse], self.ioncontdata['Pseudo'][p][:npse]) * (Coutarraynl[p])

    self.ebins_checksum = self.ebins_checksum
    self.spectrum = spec
    self.spectrum_ready = True
    if self.recombtype == FULL_NORMALIZED_RECOMBINATION:
      self.spectrum /= cross_section['sigma_total']
    return self.spectrum * collvelocity


  def calc_linelist(self, collenergy, collvelocity, specrange, specunit='A'):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)




    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """
    # convert input wavelength or energy range to angstroms
    wave = pyatomdb.util.convert_spec(specrange, specunit, 'A')

    cross_section = self.calc_crosssection(collenergy)
    Coutarraynl = cross_section['sigma']

    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    #create array to return spectrum
    l = numpy.zeros(len(self.ionlinedata), dtype=llistdtype)

    if len(self.ionlinedata) > 0:

      l['EPSILON'] = new_epsilon * collvelocity
      l['LAMBDA'] =  self.ionlinedata['LAMBDA']
      l['ENERGY'] = pyatomdb.const.HC_IN_KEV_A/l['LAMBDA']
      l['ELEMENT'] = self.ionlinedata['ELEMENT']
      l['ION'] = self.ionlinedata['ION']
      l['ION_DRV'] = self.ionlinedata['ION_DRV']
      l['UPPERLEV'] = self.ionlinedata['UPPERLEV']
      l['LOWERLEV'] = self.ionlinedata['LOWERLEV']

    self.linelist = l[(l['LAMBDA']>=wave[0]) & (l['LAMBDA']<=wave[1])]

    return self.linelist

  def calc_line_emissivity(self, collenergy, collvelocity, up, lo):

    Coutarraynl = numpy.zeros(len(self.ioncontdata), dtype = float)
    for sig in self.crosssectiondata:
      Cout = loginterp(collenergy, sig['E'], sig['C'])
      for b,item in enumerate(self.n):
        if self.n[b] == sig['n'] and self.l[b] == sig['l'] and self.s[b] == sig['S2p1']:
          Coutarraynl[b] = Cout

    #for every transition, multiplies the capture_frac with the emission, and sums the value. total emissivity
#    new_epsilon = numpy.sum(Coutarraynl * self.ionlinedata['Epsilon'], 1)
    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    i = numpy.where((self.ionlinedata['UpperLev'] == up) &
                    (self.ionlinedata['LowerLev'] == lo))[0]
    if len(i)==0:
      out_wv = 0.0
      emissivity = 0.0
    else:
      out_wv = self.ionlinedata['Lambda'][i[0]]
      emissivity = new_epsilon[i[0]]
    emissivity *= collvelocity

    ret = {'Lambda': out_wv,
           'Epsilon': emissivity,
           'up': up,
           'lo': lo}
    return(ret)


  def set_ebins(self, ebins, ebins_checksum=False):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """

    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()
    if ebins_checksum != self.ebins_checksum:
      self.ebins = ebins
      self.ebins_checksum = ebins_checksum

    for n in self.spectra.keys():
      for l in self.spectra[n].keys():
        for s in self.spectra[n][l].keys():

          self.spectra[n][l][s].set_ebins(ebins, ebins_checksum= ebins_checksum)


  def calc_crosssection(self, collenergy):

    #create array that goes through the sigma files for the known distribution of n with l from above
    Coutarraynl = numpy.zeros(len(self.ioncontdata), dtype = float)
    for sig in self.crosssectiondata:
      Cout = loginterp(collenergy, sig['E'], sig['C'])
      for b,item in enumerate(self.n):
        if self.n[b] == sig['n'] and self.l[b] == sig['l'] and \
           self.s[b]==sig['s2p1']:
          Coutarraynl[b] = Cout
    totalCout = sum(Coutarraynl)
    ret = {}
    ret['n'] = self.n
    ret['l'] = self.l
    ret['S'] = numpy.zeros(len(self.n), dtype=int)
    ret['S'][:] = -1
    ret['sigma'] = Coutarraynl
    ret['sigma_total'] = totalCout

    return(ret)



class CXIonSpectrum_N(CXIonSpectrum):
  """
  This is a class for n resolved Kronos data
  """
  import pyatomdb, numpy, os, hashlib

  """
  This is a class for nls resolved Kronos data
  """
  import pyatomdb, numpy, os, hashlib


  def __init__(self, Z,z1,crosssectiondata, \
               linedata,\
               contdata,\
               acxmodel = 8,\
               donor = False,\
               receivermass = False,\
               donormass = False):
    """
    initialize with some important variaset_collisionparambles

    Z : int
      recombining ion nuclear charge
    z1 : int
      recombining ion charge +1
    crosssectiondata: hdulist
      The  HDU with the cross section data
    linedata : hdulist
      The line emissivity data for this ion
    contdata : hdulist
      The continuum emissivity data for this ion
    acxmodel : int, optional (default = 8)
      The ACX model number in use.
    donor : string
      The donor symbol (e.g. H, He, or H2O)
    receivermass : float
      Mass in a.m.u. of the recombining particle
    donormass : float
      Mass in a.m.u. of the donor particle
    """

    # store all the info

    self.Z = Z
    self.z1 = z1

    self.Ip_donor = linedata[1].header['Ip_D']
    self.crosssectiondata = crosssectiondata
    self.donor = donor
    self.receivermass = receivermass
    self.donormass = donormass
    self.acxmodel = acxmodel

    # set up some empty variables, to be fleshed out
    self.ebins_checksum = "0" # store hash of ebins
    # find, and store, the relevant HDUs from the line and continuum datafiles.

    #find the hdu location for the matching Z and z1 for line and continuum

    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0][0]

    self.ionlinedata = linedata[i+2].data

    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0][0]

    self.ioncontdata = contdata[i+2].data

    #find the n l s values from the header and write to arrays
    self.n = numpy.zeros(len(self.ioncontdata))
    self.l = numpy.zeros(len(self.ioncontdata))
    self.s = numpy.zeros(len(self.ioncontdata))

    for j,item in enumerate(self.ioncontdata):
       headstr = "CONT_*"
       headstr = headstr.replace("*", str(j))
       nlsstr = contdata[i+2].header[headstr]
       nlslist = nlsstr.replace('n',' ').replace('l',' ').replace('s',' ').split()

       self.n[j] = nlslist[0]
       self.l[j] = nlslist[1]
       self.s[j] = nlslist[2]

    #to store the location of all the HDU data
    self.linedata = linedata
    self.contdata = contdata

    self.spectra ={}
    self.ebins=0.0
    self.ebins_checksum="0"



  def calc_spectrum(self, collenergy, collvelocity, Tbroaden, vbroaden, ebroaden):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    linedata : hdulist
      The line emissivity data
    contdata : hdulist
      The continuum emissivity data
    acxmodel : int
      n, l shell distribution, between 1 and 8


    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """


    # This will be returned NOT accounting for ionization balance or
    # elemental abundance. It is assumed that this will be handled at
    # a higher level, and indeed, since this call is quite time consuming,
    # it shouldn't be made at all for unnecessary ions.


    # cycle through every single entry in cross section file
    #sigmadata = self.crosssectiondata

    #find energy range
    #emin = min(sigmadata['E_min'])
    #emax = min(sigmadata['E_max'])
#    print("N Spectrum calc")

    cross_section = self.calc_crosssection(collenergy)

    Coutarraynl = cross_section['sigma']

    #for each shell, multiplies the capture_frac with the emission, and sums the value. total emissivity
#    new_epsilon = numpy.sum(Coutarraynl * self.ionlinedata['Epsilon'], 1)
    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)
#    new_epsilon = numpy.zeros(len(self.ionlinedata))
#    for iline,item in enumerate(self.ionlinedata):
#      new_epsilon[iline] = numpy.sum(Coutarraynl * self.ionlinedata['Epsilon'][iline,:])

    #set up array to return spectrum
    spec = numpy.zeros(len(self.ebins)-1, dtype=float)

    #######line broadening
    self.broaden_limit = 1e-20
    self.broaden_object = _Gaussian_CDF()

    if ((Tbroaden <= 0.0) & \
          (vbroaden <=0.0)):

      if len(self.ionlinedata) > 0:
        #turn line wavelengths into energies
        en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
        #puts lines into the right spectral bins using epsilon array
        a,b= numpy.histogram(en, bins = self.ebins, weights = new_epsilon)
        #add the line emission per bin to the spectrum
        spec += a

    elif len(self.ionlinedata) > 0:

      en = pyatomdb.const.HC_IN_KEV_A/self.ionlinedata['Lambda']
      ind = new_epsilon>self.broaden_limit
      nonind = ~ind

      # calculate the widths of the strong lines
      llist = self.ionlinedata[ind]

      # get a raw dictionary of masses in amu
      masslist = pyatomdb.atomic.Z_to_mass(1,raw=True)

      #assign values of temperature and velocity broadening
      T = Tbroaden
      Tb = pyatomdb.util.convert_temp(T, 'K','keV')*pyatomdb.const.ERG_KEV/(masslist[llist['Element']]*1e3*pyatomdb.const.AMUKG)

      if vbroaden <0:
        vbroaden = 0.0
        vb=0.0
      else:
        vb = (vbroaden * 1e5)**2

      wcoeff = numpy.sqrt(Tb+vb) / (pyatomdb.const.LIGHTSPEED*1e2)

      elines = en[ind]
      n_e = new_epsilon[ind]
      width = wcoeff*elines

      # Filter out lines more than NSIGMALIMIT sigma outside the range
      NSIGMALIMIT=4
      eplu = elines+NSIGMALIMIT*width
      eneg = elines-NSIGMALIMIT*width
      emax = max(self.ebins)
      emin = min(self.ebins)
      # identify all the good lines!
      igood = numpy.where(((elines >= emin) & (eneg < emax))  |\
                ((elines < emax) & (eplu > emin)))[0]

      spec = numpy.zeros(len(self.ebins))

      for iline in igood:
        #if width
        spec += self.broaden_object.broaden(elines[iline],\
                         width[iline],self.ebins)*n_e[iline]
      spec=spec[1:]-spec[:-1]
      s,z = numpy.histogram(en[nonind], \
                            bins = self.ebins,\
                            weights = new_epsilon[nonind])
      spec+=s


    #######

    #add to spectrum for the Coutarraynl n l s values
    if len(self.ioncontdata) > 0:
      for p,item in enumerate(self.ioncontdata):
        if Coutarraynl[p] > 0.0:
          if self.ioncontdata['N_cont'][p] > 2:
            ncont = self.ioncontdata['N_Cont'][p] #shorten array
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Cont'][p][:ncont], self.ioncontdata['Continuum'][p][:ncont]) * (Coutarraynl[p])

          if self.ioncontdata['N_Pseudo'][p] > 2:
            npse = self.ioncontdata['N_Pseudo'][p]
            spec += self.expand_E_grid(self.ebins, self.ioncontdata['E_Pseudo'][p][:npse], self.ioncontdata['Pseudo'][p][:npse]) * (Coutarraynl[p])

    self.ebins_checksum = self.ebins_checksum
    self.spectrum = spec
    self.spectrum_ready = True

    # renormalize by total cross section
    if self.recombtype == FULL_NORMALIZED_RECOMBINATION:
      if cross_section['sigma_total'] > 0:
        print("rescaling spectrum N by %e"%(cross_section['sigma_total']))

        self.spectrum/=cross_section['sigma_total']
    return self.spectrum * collvelocity


  def calc_linelist(self, collenergy, collvelocity, specrange, specunit='A'):
    """
    Calculate the spectrum of the data

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins (in keV) to calcualte the spectrum on
    collenergy : flaot
      The collision energy (keV/amu)
    collvelocity : float
      The velocity of the center of mass (cm/s)
    specrange : [float, float]
      Minimum and maximum values for interval in which to search
    specunit : {'Angstrom','keV'}
      Units for specrange (default A)




    RETURNS
    -------
    emissivity : array(float)
      Emissivity * velocity in photons cm4 s-2

    """
    # convert input wavelength or energy range to angstroms
    wave = pyatomdb.util.convert_spec(specrange, specunit, 'A')

    cross_section = self.calc_crosssection(collenergy)
    Coutarraynl = cross_section['sigma']

    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)

    #create array to return spectrum
    l = numpy.zeros(len(self.ionlinedata), dtype=llistdtype)

    if len(self.ionlinedata) > 0:

      l['EPSILON'] = new_epsilon * collvelocity
      l['LAMBDA'] =  self.ionlinedata['LAMBDA']
      l['ENERGY'] = pyatomdb.const.HC_IN_KEV_A/l['LAMBDA']
      l['ELEMENT'] = self.ionlinedata['ELEMENT']
      l['ION'] = self.ionlinedata['ION']
      l['ION_DRV'] = self.ionlinedata['ION_DRV']
      l['UPPERLEV'] = self.ionlinedata['UPPERLEV']
      l['LOWERLEV'] = self.ionlinedata['LOWERLEV']

    self.linelist = l[(l['LAMBDA']>=wave[0]) & (l['LAMBDA']<=wave[1])]

    return self.linelist


  def calc_line_emissivity(self, collenergy, collvelocity, up, lo):

    l = -1 * (self.acxmodel%4)
    # if l = -4, this returns 0. Correct
    if l==0: l=-4

    #create array that goes through the sigma files for the known distribution of n with l from above
    Coutarraynl = numpy.zeros(len(self.ioncontdata), dtype = float)
    for sig in self.crosssectiondata:
      Cout = loginterp(collenergy, sig['E'], sig['C'])
      for b,item in enumerate(self.n):
        if self.n[b] == sig['n'] and self.l[b] == l:
          Coutarraynl[b] = Cout

    #for each shell, multiplies the capture_frac with the emission, and sums the value. total emissivity
#    new_epsilon = numpy.sum(Coutarraynl * self.ionlinedata['Epsilon'], 1)
    if len(Coutarraynl)==0:
      new_epsilon = numpy.zeros(len(self.ionlinedata))

    else:
      new_epsilon = numpy.sum(self.ionlinedata['Epsilon']*Coutarraynl,1)


    i = numpy.where((self.ionlinedata['UpperLev'] == up) &
                    (self.ionlinedata['LowerLev'] == lo))[0]
    if len(i)==0:
      out_wv = 0.0
      emissivity = 0.0
    else:
      out_wv = self.ionlinedata['Lambda'][i[0]]
      emissivity = new_epsilon[i[0]]

    emissivity *= collvelocity
    ret = {'Lambda': out_wv,
           'Epsilon': emissivity,
           'up': up,
           'lo': lo}

    return ret



  def set_ebins(self, ebins, ebins_checksum=False):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """

    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()
    if ebins_checksum != self.ebins_checksum:
      self.ebins = ebins
      self.ebins_checksum = ebins_checksum

    for n in self.spectra.keys():
      for l in self.spectra[n].keys():
        self.spectra[n][l].set_ebins(ebins, ebins_checksum= ebins_checksum)

  def calc_crosssection(self, collenergy):

    # there is no l distribution, so use ACX1 settings
    l = -1 * (self.acxmodel%4)
    # if l = -4, this returns 0. Correct back to -4
    if l==0: l=-4

    #create array that goes through the sigma files for the known distribution of n with l from above
    Coutarraynl = numpy.zeros(len(self.ioncontdata), dtype = float)
    for sig in self.crosssectiondata:
      Cout = loginterp(collenergy, sig['E'], sig['C'])
      for b,item in enumerate(self.n):
        if self.n[b] == sig['n'] and self.l[b] == l:
          Coutarraynl[b] = Cout
    totalCout = sum(Coutarraynl)
    ret = {}
    ret['n'] = self.n
    ret['l'] = self.l
    ret['S'] = numpy.zeros(len(self.n), dtype=int)
    ret['S'][:] = -1
    ret['sigma'] = Coutarraynl
    ret['sigma_total'] = totalCout

    return(ret)


class CXShellSpectrum():
  """
  Holds a single n, l, shell spectrum
  """
  import pyatomdb, numpy, os, hashlib

  import scipy,scipy.integrate

  def __init__(self, Z, z1, n, l, \
               linedata,\
               contdata, s = False):
    """
    Initialization Function

    PARAMETERS
    ----------
    Z : intspitzencandidate
      atomic number of recombining ion
    z1 : int
      ion charge +1 of recombining ion
    n : int
      n shell of this spectrum
    l : int
      l shell of this spectrum (-ve numbers denote distributions)
    linedata : HDU
      The HDU for the line emissivity for recombination into this shell
    contdata : HDU
      The HDU for the continuum emissivity due acx2_xspec.py", line 1835, in calc_spectrum
    en = pyatomdb.const.HC_IN_KEV_A/self.linedata['Lambda']
to recombination into this shell.
    s : int, optional
      The 2S+1 of the quantum number

    RETURNS
    -------
    None
    """

    self.Z = Z
    self.z1 = z1
    self.n = n
    self.l = l
    self.s = s
    self.linedata = linedata
    self.contdata = contdata
    self.ebins_checksum = "0"
    self.spectrum_ready=False
    self.spectrum = 0.0



  #-------------------------------------------------------------------------------
  def expand_E_grid(self, eedges, Econt_in_full, cont_in_full):

    """
    Code to expand the compressed continuum onto a series of bins.

    Parameters
    ----------
    eedges : float(array)
      The bin edges for the spectrum to be calculated on, in units of keV
    Econt_in_full: float(array)
      The compressed continuum energies (keV)
    cont_in_full: float(array)
      The compressed continuum emissivities (ph cm3 s-1 keV-1)

    Returns
    -------
    float(array)
      len(bins)-1 array of continuum emission, in units of \
      photons cm^3 s^-1 bin^-1
    """

  #  History
  #  -------
  #  Version 0.1 - initial release
  #    Adam Foster July 17th 2015

    import scipy.integrate
    cont_in = cont_in_full
    Econt_in = Econt_in_full
    n = len(cont_in_full)

    # Append desired bin edges to the end of this array
    E_all = numpy.append(Econt_in, eedges)

    # interpolate the continuum at these bin edges
    cont_tmp = numpy.interp(eedges, Econt_in, cont_in)

    # add continuum at bin edges to the end
    C_all = numpy.append(cont_in, cont_tmp)

    # get sorted indexes for energy
    iord = numpy.argsort(E_all)

    # reorder the arrays in energy order
    E_all = E_all[iord]
    C_all = C_all[iord]

    # pull out the points which are from the output grid (iord > npoints)
    ihi = numpy.where(iord>=n)[0]

    # do cumulative sum integration
    cum_cont = scipy.integrate.cumulative_trapezoid(C_all, E_all, initial=0)

    # extract the cumulative sum at the output points
    C_out = cum_cont[ihi]

    # convert to flux per bin
    cont = C_out[1:]-C_out[:-1]

    return cont





  def calc_spectrum(self, ebins, ebins_checksum, Tbroaden, vbroaden, ebroaden):
    import scipy,scipy.integrate

    if self.ebins_checksum == ebins_checksum:
      return self.spectrum
    else:
      # line emission
      spec = numpy.zeros(len(ebins)-1, dtype=float)

      if len(self.linedata) > 0:

        en = pyatomdb.const.HC_IN_KEV_A/self.linedata['Lambda']

        a,b= numpy.histogram(en, bins = ebins, weights = self.linedata['Epsilon'])

        spec += a
      if len(self.contdata) > 0:
        if self.contdata['N_cont'][0] > 2:

          spec += self.expand_E_grid(ebins, self.contdata['E_Cont'][0], self.contdata['Continuum'][0])
        if self.contdata['N_Pseudo'][0] > 2:
          spec += self.expand_E_grid(ebins, self.contdata['E_Pseudo'][0], self.contdata['Pseudo'][0])

      self.ebins_checksum = ebins_checksum
      self.spectrum = spec
      self.spectrum_ready = True

      return self.spectrum




  def calc_line_emissivity(self, up, lo):
    import scipy,scipy.integrate

    eps=0.0
    wv=0.0

    if len(self.linedata) > 0:

      en = pyatomdb.const.HC_IN_KEV_A/self.linedata['Lambda']

      l = self.linedata[ (self.linedata['UpperLev']==up) & (self.linedata['LowerLev']==lo) ]

      for ll in l:
        wv = ll['Lambda']
        eps += ll['Epsilon']

    self.line_emissivity={'Lambda':wv, 'Epsilon':eps, 'up':up, 'lo':lo}

    return self.line_emissivity


  def set_ebins(self, ebins, ebins_checksum):
    """
    Set the energy bins, also in each donor model

    PARAMETERS
    ----------
    ebins : array(float)
      The energy bins in keV for the spectrum

    RETURNS
    -------
    None
    """

    if ebins_checksum == False:
      ebins_checksum = hashlib.md5(ebins).hexdigest()

    if self.ebins_checksum != ebins_checksum:
      # remove the spectrum
      self.spectrum=0.0
      self.spectrum_ready=False




class DummyCXShellSpectrum():
  """
  Placeholder for a blank spectrum
  """
  def __init__(self, Z, z1, n, l, s=False):
    """
    Initialization Function

    PARAMETERS
    ----------
    Z : intspitzencandidate
      atomic number of recombining ion
    z1 : int
      ion charge +1 of recombining ion
    n : int
      n shell of this spectrum
    l : int
      l shell of this spectrum (-ve numbers denote distributions)

    RETURNS
    -------
    None
    """

    self.Z = Z
    self.z1 = z1
    self.n = n
    self.l = l
    self.s = s


  def calc_spectrum(self, ebins, ebins_checksum, Tbroaden, vbroaden, ebroaden):
    """
    Return zeros

    PARAMETERS
    ----------
    ebinshash : string
      md5 sum of ebins
    ebins : array(float)
      array of floats

    RETURNS
    -------
    0.0 : there is no data here, so return zero.
    """

    return 0.0

  def set_ebins(self, ebins, ebins_checksum=False):
    """
    Dummy

    PARAMETERS
    ----------
    ebinshash : string
      md5 sum of ebins
    ebins : array(float)
      array of floats

    RETURNS
    -------
    None
    """

    pass



class _Gaussian_CDF():
  """
  For fast interpolation, pre-calculate the CDF and interpolate it when
  broadening lines

  Parameters
  ----------
  None

  Examples
  --------

  Create a CDF instance:

  >>> s=_Gaussian_CDF()

  Broaden a line on ebins grid, with centroid and width.

  >>> cdf = a.broaden(centroid, width, ebins)

  Convert to flux in each bin

  >>> flux= cdf[1:]-cdf[:-1]


  """
  def __init__(self):
    from scipy.stats import norm
    self.x = numpy.linspace(-6,6,2400)
    self.cdf = norm.cdf(self.x)
    self.broadentype='Gaussian'

  def broaden(self, centroid, width, ebins):
    """
    Broaden a line, return CDF

    Parameters
    ----------
    centroid : float
      The line energy (keV)
    width : float
      The sigma of the normal distribution, in keV
    ebins : array(float)
      Energy grid to return CDF on

    Returns
    -------
    CDF : array(float)
      cumulative flux distribution of linen at each bin edge.
    """

    # move the energy grid
    etmp = (ebins-centroid)/width

    # interpolate to get the appropriate CDF values
    ret=numpy.interp(etmp, self.x, self.cdf)

    return ret


def lvdat_to_termsymb(d):
  LLIST = 'SPDFGHIJKLMNOQRTUVWZYZabcdefghijklmnopqrstuvwxyz'

  s2p1 = int(d['S_QUAN']*2 +1)
  ssym = "%i"%(s2p1)
  lquan = d['L_QUAN']


  if lquan ==-1:
    lsym = '-1'
  else:
    lsym = LLIST[lquan]

  j2 = d['LEV_DEG']-1
  if j2%2 == 0:
    jsym = '%i'%(j2/2)
  else:
    jsym = '%i/2'%(j2)
  dout = "%s(%s)%s"%(ssym, lsym, jsym)

  return(dout)

def get_full_line_info(linelist, filemap="$ATOMDB/filemap_acx", atomdbroot="$ATOMDB"):
  """
  Gets the data for the full linelist
  """


  datacache = {}
  settings = {}
  settings['filemap'] = filemap
  settings['atomdbroot']=atomdbroot

  dout = numpy.zeros(len(linelist), dtype = llistfulldtype)

  for k in ['EPSILON','LAMBDA','ENERGY',\
                                   'ELEMENT','ION','ION_DRV',\
                                   'UPPERLEV','LOWERLEV']:
    dout[k] =linelist[k]

  for i, l in enumerate(linelist):
    lvdat = pyatomdb.atomdb.get_data(l['ELEMENT'], l['ION'], 'LV', datacache=datacache, \
             settings=settings)

    ladat = pyatomdb.atomdb.get_data(l['ELEMENT'], l['ION'], 'LA', datacache=datacache, \
             settings=settings)

    dout[i]['UPPERCFG'] = lvdat[1].data['ELEC_CONFIG'][l['UPPERLEV']-1]
    dout[i]['UPPERSYM'] = lvdat_to_termsymb(lvdat[1].data[l['UPPERLEV']-1])

    dout[i]['LOWERCFG'] = lvdat[1].data['ELEC_CONFIG'][l['LOWERLEV']-1]
    dout[i]['LOWERSYM'] = lvdat_to_termsymb(lvdat[1].data[l['LOWERLEV']-1])

    j = numpy.where( (ladat[1].data['UPPER_LEV']==l['UPPERLEV']) &\
                     (ladat[1].data['LOWER_LEV']==l['LOWERLEV']))[0][0]
    dout[i]['EINSTEIN_A'] = ladat[1].data['EINSTEIN_A'][j]
    dout[i]['B_RATIO'] = ladat[1].data['EINSTEIN_A'][j]/lvdat[1].data['ARAD_TOT'][l['UPPERLEV']-1]

  return(dout)


def test():
  # this is a test routine. Let's see what we can do!

  myacx = ACXModel()

  myacx.add_donor('H', \
                'tmp_cxline2.fits', \
                'tmp_cxcont.fits', \
                'tmp_cx_sigma.fits',\
                elements=[1,2,8])
    ### TESTFUDGE
#  zzz=input('woo2')

  myacx.set_temperature(5.0, teunits='keV')



 # myacx.velocity = 500*1000*100. # in cm/s
  elo = 12.398425/310
  ehi = 12.398425/10
  ebins = numpy.linspace(elo, ehi, 1000)
  myacx.set_ebins(ebins)
  myacx.set_collisiontype(4, 'cm/s')
  return myacx


# def download_acx_emissivity_files(adbroot=None, userid, version=None):

  # """
  # Download the AtomDB CX model emissivity files for AtomDB"

  # This code will go to the AtomDB FTP site and download the necessary files.
  # It will then unpack them into a directory adbroot. It will not
  # overwrite existing files with the same md5sum (to avoid pointless updates)
  # but it will not know this until it has downloaded and unzipped the main
  # file.

  # Parameters
  # ----------

  # adbroot : string
    # The location to install the data. Typically should match $ATOMDB
  # userid : string
    # An 8 digit ID number. Usually passed as a string, but integer
    # is also fine (provided it is all numbers)
  # version : string
    # The version string, including donors for the release, e.g. "2_1_0-h_he"

  # Returns
  # -------
  # None
  # """

  # import tarfile

  # # set up remote file name
  # fname = "acx2-v%s-data.tar.bz2"%(version)

  # if adbroot is None:
    # adbroot = os.environ.get('ATOMDB')
    # if adbroot is None:
      # print("Error: to automatically download and install AtomDB files, must have $ATOMDB environment variable set")

  # a=curl.Curl()
  # version=a.get('%s/releases/LATEST_ACX2'%(const.FTPPATH))[:-1].decode(encoding='ascii')
  # a.close()

  # # set up temporary directory to hold data

  # if adbroot[0] != '/':
    # # is a relative path
    # adbroot = "%s/%s"%(os.getcwd(), adbroot)

  # mkdir_p(adbroot)
  # if adbroot[-1]=='/':
    # tmpdir = adbroot+'installtmp'
  # else:
    # tmpdir = adbroot+'/installtmp'


  # print("making directory %s"%(tmpdir))
  # mkdir_p(tmpdir)

  # # get the files
  # urllib.request.urlcleanup()
  # fnameout = wget.download('%s/releases/%s'%(pyatomdb.const.FTPPATH,fname), out="%s/%s"%(tmpdir, fname))
  # # collect user statistics if allowed.
  # record_upload(fname)

  # #uncompress
  # print("")
  # print("Uncompressing ... ", end=' ')

  # tf = tarfile.open(name=fnameout, mode='r:bz2')
  # tf.extractall(path=tmpdir)
  # print("Uncompressed")
  # # copy the files
  # dirname = 'acx2_v%s-data'%(version)
  # for l in os.listdir('%s/%s'%(tmpdir, dirname)):
    # print("moving %s/%s/%s to %s/%s"%(tmpdir, dirname, l, adbroot, l))
    # shutil.move("%s/%s/%s"%(tmpdir, dirname, l), "%s/%s"%(adbroot, l))

  # print("...done")

  # shutil.rmtree(tmpdir)



