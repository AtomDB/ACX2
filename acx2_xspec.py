import pyatomdb, numpy, os, hashlib

# this is the library for running the ACX2 model.

#some constants!

CENTEROFMASS=1
UNIVERSAL_CX_CROSSSECTION=3e-15 # in cm^2
DEBUG=False

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

  def __init__(self):
    """
    Initialize the ACX Model
    """

    self.DonorList = []
    self.ebins = False
    self.ebins_checksum = False
    self.temperature = False
    self.ebins_set = False
    self.temperature_set = False
    self.collision_set = False
    self.abund = {}
    for Z in range(1,31):
      self.abund[Z] = 1.0


  def add_donor(self, donor, \
                donor_linefile, \
                donor_contfile, \
                donor_crosssectionfile,\
                abundset='AG89',\
                elements=list(range(1,31))):
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

    self.DonorList.append(ACXDonorModel(donor, donor_linefile,\
                                        donor_contfile,\
                                        donor_crosssectionfile,\
                                        abundset=abundset,\
                                        elements=elements))
    if self.ebins_set:
      self.DonorList[-1].set_ebins(ebins)

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

    for donor in self.DonorList:
      try:
        i = donorlist.index(donor.lower())
        donor.set_donorabund(abundlist[i])
      except ValueError:
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
    collunits : string {'cm/s', 'kev/amu','kev/u', 'ev/amu','ev/u'}
      The units of collvalue (case insensitive)
    """
#    print("Collunits", collunits)
    self.collisionunits=collunits
    self.collisiontype=colltype
    self.collisionset=True
    for donor in self.DonorList:
      donor.set_collisiontype(colltype, collunits=collunits)








  def calc_spectrum(self, collvalue):
    """
    Calculate the spectrum for all the donors, sum.

    PARAMETERS
    ----------
    collvalue : float
      The collision parameter (kev/u or cm/s, depending) to calculate the spectrum

    RETURNS
    -------
    ret : array(float)
      The spectrum, in ph cm^3 s-1 bin-1
    """
    retset=False

    if DEBUG:
      ret = []
      for donor in self.DonorList:
        ret.append(donor.calc_spectrum(collvalue))

      return ret
    else:
      ret=False
      for donor in self.DonorList:

        if retset==True:
          ret += donor.calc_spectrum(collvalue)
        else:
          ret = donor.calc_spectrum(collvalue)
          retset=True

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
        print("Error: specified incompatible number of elements from set_abund")
        return

    for i in range(len(abundvec)):
      self.abund[elements[i]] = abundvec[i]

    for donor in self.DonorList:
      donor.set_abund(abundvec, elements=elements)



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
               donor_linefile, \
               donor_contfile, \
               donor_crosssectionfile,\
               abundset='AG89',\
               elements=list(range(1,31)),\
               acxmodel = 8):

    self.donor = donor.lower() # store in lower case
    self.donor_linefile = os.path.expandvars(donor_linefile)
    self.donor_contfile = os.path.expandvars(donor_contfile)
    self.donor_crosssectionfile = os.path.expandvars(donor_crosssectionfile)
#    input('DON1')
    try:
      self.linedata = pyatomdb.pyfits.open(self.donor_linefile)
    except:
      print("Cannot open line data file %s"%(self.linedata))
      raise

    try:
      self.contdata = pyatomdb.pyfits.open(self.donor_contfile)
    except:
      print("Cannot open continuum data file %s"%(self.contdata))
      raise

    try:
      self.crosssectiondata = pyatomdb.pyfits.open(self.donor_crosssectionfile)
      self.donormass = self.crosssectiondata[1].header['DonMass']

    except:
      print("Cannot open cross section data file %s"%(self.crosssectiondata))
      raise

    # create a structure for the spectral data
    self.spectra = {}

    # set the default abundance of the donor to 1.0
    self.set_donorabund(1.0)


    # create an abundance vecotr

    self.abundset=abundset
    self.default_abundset=abundset
#    input('DON2')


    # if elements are specified, use them. Otherwise, use Z=1-30
    if pyatomdb.util.keyword_check(elements):
      self.elements = elements
    else:
      self.elements=list(range(1,31))

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

    # temporary things
    self.datacache={}

    self.ebins = 0.0
    self.ebins_checksum = ''


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
#    input("DON2.2.1")
    # read in the abundance the raw data was calculated on
    old = pyatomdb.atomdb.get_abundance(abundset=self.default_abundset)
    #input("old = ", old)
#    input("DON2.2.2")

    # read in the new abundance
    new = pyatomdb.atomdb.get_abundance(abundset=abundstring)
#    input("DON2.2.3")

    # divide the 2, store the replacement ratio to self.abundsetvector
    for Z in list(self.abundsetvector.keys()):
      self.abundsetvector[Z]=new[Z]/old[Z]
#    input("DON2.2.4")

    # update the current abundance string to represent your input
    self.abundset=abundstring
#    input("DON2.2.5")

  #def set_abundance(self, abundance):loginterp


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

#      self.crosssectiondata = crosssectiondata[ihdu+2].data
#      self.coupling = crosssectiondata['INDEX'].data['resn'][i].decode('ascii')
#      self.DonorMass = crosssectiondata[ihdu+2].header['DONMASS']
#      self.RecvMass = crosssectiondata[ihdu+2].header['RECMASS']

    else:
      ihdu = -1
      resolution = 'ACX1'
#      self.RecvMass=pyatomdb.atomic.Z_to_mass(self.Z)
#      self.DonorMass=crosssectiondata['INDEX'].header['DONMASS']
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
          self.collvelocity[Z] = 1e5* numpy.sqrt(4786031.3*self.collenergy[Z]*1e-3/25.)
        else:
          print("set_collisionparam: unknown collision units %s"%(self.collisionunits))
          return
      elif self.collisiontype==2:
        if self.collisionunits.lower() == 'cm/s':
          # this is the reduced mass velocity
          self.collvelocity[Z] = collisionparam*1.0
          self.collenergy[Z] = 25 * (collisionparam/1e5)**2/4786031.3
      elif self.collisiontype == 3:
          # donor ion is moving
        if self.collisionunits.lower() == 'cm/s':
          self.collvelocity[Z] = self.donormass*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
          self.collenergy[Z] = 25 * (self.collvelocity[Z]/1e5)**2/4786031.3
      elif self.collisiontype == 4:
          # receiver ion is moving
        if self.collisionunits.lower() == 'cm/s':
          self.collvelocity[Z] = pyatomdb.atomic.Z_to_mass(Z)*collisionparam/(self.donormass+pyatomdb.atomic.Z_to_mass(Z))
          self.collenergy[Z] = 25 * (self.collvelocity[Z]/1e5)**2/4786031.3

      else:
        print("*** ERROR: Unknown collision unit %s: should be kev/amu or cm/s ***" %(self.collisionunits))
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
#          input('ACXDM_fuckit1[%i][%i]'%(Z,z1))

          self.spectra[Z][z1].set_ebins(self.ebins, ebins_checksum=self.ebins_checksum)


#  def calc_spectrum(self, ebins, ionfrac, collparam, linedata, contdata, ebins_checksum=False):
  def calc_spectrum(self, collparam):
    """
    Calculate the spectrum we want

    PARAMETERS
    ----------

    ebins : array(float)
      Energy bin edges in keV
    ionfrac : dict
      Ionization balance
    collparam : float
      Collision energy, or velocity, as determined by set_collisiontype, in kev/u or cm/s
    linedata : HDUlist
      The CX line data
    contdata : HDUlist
      The CX continuum data


    """

    # get the velocity etc
    self.set_collisionparam(collparam)


    # set up return array for spectrum
    if DEBUG:
      self.emiss_debug = {}
    self.emiss = numpy.zeros(len(self.ebins)-1, dtype=float)


    for Z in self.elements:
      if self.abund[Z] > 0.0:
        if not Z in self.spectra.keys():
          self.spectra[Z]={}
        for z1 in range(2, Z+2):# z1 here is the recombin*ing* ion charge +1
          if self.abund[Z]*self.ionfrac[Z][z1-1] > 1e-10: #ionfrac is indexed from 0, not 1
#            print("Z=%i, z1=%i"%(Z, z1))


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

            self.spectra[Z][z1].set_ebins(self.ebins, ebins_checksum=self.ebins_checksum)

            if DEBUG:
              if not Z in self.emiss_debug.keys():
                self.emiss_debug[Z] = {}

#              self.emiss_debug[Z][z1] = self.spectra[Z][z1].calc_spectrum(self.ebins, self.collenergy[Z], self.collvelocity[Z], self.linedata, self.contdata, self.acxmodel)
              self.emiss_debug[Z][z1] = self.spectra[Z][z1].calc_spectrum(self.collenergy[Z], self.collvelocity[Z])

              self.emiss +=  self.emiss_debug[Z][z1] * self.abund[Z] * self.ionfrac[Z][z1-1]

            else:

#              self.emiss += self.spectra[Z][z1].calc_spectrum(self.ebins, self.collenergy[Z], self.collvelocity[Z], self.linedata, self.contdata, self.acxmodel) *\
#                     self.abund[Z] * self.ionfrac[Z][z1-1]
#              print("Z",Z,"z1",z1)
#              print("self.abund",self.abund)
#              print("self.ionfrac",self.ionfrac)
#              print("self.spectra",self.spectra.keys())
#              print("self.spectra[%i][%i]"%(Z,z1),self.spectra[Z][z1])
#              print("self.collenergy[Z]",self.collenergy[Z])
#              print("self.collvelocity[Z]",self.collvelocity[Z])
              self.emiss += self.spectra[Z][z1].calc_spectrum(self.collenergy[Z], self.collvelocity[Z]) *\
                     self.abund[Z] * self.ionfrac[Z][z1-1]


    return self.emiss


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
          self.ionfrac[Z] = pyatomdb.atomdb.apec.solve_ionbal_eigen(Z, self.temperature,\
                                                                   teunit='kev', datacache=self.datacache)
    self.ionfrac_from_temperature = True


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

  def set_abund(self, abund, elements=False):
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
        print("Error: specified incompatible number of elements from set_abund")
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
      self.RecvMass=pyatomdb.atomic.Z_set_collisionparamto_mass(self.Z)
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

    self.linedataindexes={} # to store the location of all the HDU data


    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0]


    nlist = pyatomdb.util.unique(linedata[1].data['n'][i])
    nlist.sort()
    for n in nlist:
       self.linedataindexes[n] = {}

       ii = numpy.where((linedata[1].data['Z'] == self.Z) &\
                        (linedata[1].data['z1'] == self.z1) &\
                        (linedata[1].data['n'] == n))[0]
       if len(ii) > 0:
         llist = pyatomdb.util.unique(linedata[1].data['l'][ii])
         for l in llist:
           self.linedataindexes[n][l] = numpy.where((linedata[1].data['Z'] == self.Z) &\
                        (linedata[1].data['z1'] == self.z1) &\
                        (linedata[1].data['n'] == n) &\
                        (linedata[1].data['l'] == l))[0][0]+2



    self.contdataindexes={} # to store the location of all the HDU data


    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0]


    nlist = pyatomdb.util.unique(contdata[1].data['n'][i])
    nlist.sort()
    for n in nlist:
       self.contdataindexes[n] = {}

       ii = numpy.where((contdata[1].data['Z'] == self.Z) &\
                        (contdata[1].data['z1'] == self.z1) &\
                        (contdata[1].data['n'] == n))[0]
       if len(ii) > 0:
         llist = pyatomdb.util.unique(contdata[1].data['l'][ii])
         for l in llist:
           self.contdataindexes[n][l] = numpy.where((contdata[1].data['Z'] == self.Z) &\
                        (contdata[1].data['z1'] == self.z1) &\
                        (contdata[1].data['n'] == n) &\
                        (contdata[1].data['l'] == l))[0][0]+2

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

#  def calc_spectrum(self, ebins, collenergy, collvelocity, linedata, contdata, acxmodel):
  def calc_spectrum(self, collenergy, collvelocity):
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
#    print("nprime = %f, n="%(nprime), n,", nfrac=",nfrac)
    # set up return arrays


#    ebinshash = hashlib.md5(ebins).hexdigest() # tracks for changes to ebins

    emissivity = numpy.zeros(len(self.ebins)-1, dtype=float)


#    print("Linedataindexes:", self.linedataindexes)

    for inn, nn in enumerate(n):
      if not nn in self.spectra.keys():
        self.spectra[nn] = {}
#      print("n=%i, l=%i"%(nn,l))
      if not l in self.spectra[nn].keys():


        try:
          self.spectra[nn][l] = CXShellSpectrum(self.Z, self.z1, nn, l, \
                                               self.linedata[self.linedataindexes[nn][l]].data,\
                                               self.contdata[self.contdataindexes[nn][l]].data)
        except KeyError:
          # Case where there is no data for this Z, z1, nn, l
          self.spectra[nn][l] = DummyCXShellSpectrum(self.Z, self.z1, nn, l)


      # do the actual calculation, add to return value
      emissivity += nfrac[inn] * self.spectra[nn][l].calc_spectrum(self.ebins, self.ebins_checksum) * UNIVERSAL_CX_CROSSSECTION

#    print("sum emissivity:", sum(emissivity))
    # need velocity in cm s-1
    return emissivity * collvelocity



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

    self.linedataindexes={} # to store the location of all the HDU data
    self.linedata=linedata # to store the location of all the HDU data


    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0]


    for ii in i:
      n = linedata[1].data['n'][ii]
      l = linedata[1].data['l'][ii]
      s = linedata[1].data['S2p1'][ii]

      if not n in self.linedataindexes.keys():
        self.linedataindexes[n] ={}
      if not l in self.linedataindexes[n].keys():
        self.linedataindexes[n][l] ={}

      self.linedataindexes[n][l][s] = ii+2


    self.contdataindexes={} # to store the location of all the HDU data
    self.contdata=contdata # to store the location of all the HDU data


    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0]


    for ii in i:
      n = contdata[1].data['n'][ii]
      l = contdata[1].data['l'][ii]
      s = contdata[1].data['S2p1'][ii]

      if not n in self.contdataindexes.keys():
        self.contdataindexes[n] ={}
      if not l in self.contdataindexes[n].keys():
        self.contdataindexes[n][l] ={}

      self.contdataindexes[n][l][s] = ii+2

    self.spectra ={}
    self.ebins=0.0
    self.ebins_checksum="0"






  def calc_spectrum(self, collenergy, collvelocity):
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

    emissivity = numpy.zeros(len(self.ebins)-1, dtype=float)
#    ebinshash = hashlib.md5(self.ebins).hexdigest() # tracks for changes to ebins


    for sig in self.crosssectiondata:
      Cout = loginterp(collenergy, sig['E'], sig['C'])

      if Cout > 0:

        if not sig['n'] in self.spectra.keys():
          self.spectra[sig['n']] = {}
        if not sig['l'] in self.spectra[sig['n']].keys():
          self.spectra[sig['n']][sig['l']] = {}
        if not sig['S2p1'] in self.spectra[sig['n']][sig['l']].keys():
          try:
            self.spectra[sig['n']][sig['l']][sig['S2p1']] = \
               CXShellSpectrum(self.Z, self.z1, sig['n'], sig['l'], \
               self.linedata[self.linedataindexes[sig['n']][sig['l']][sig['S2p1']]].data,\
               self.contdata[self.contdataindexes[sig['n']][sig['l']][sig['S2p1']]].data)
          except KeyError:
          # Case where there is no data for this Z, z1, nn, l
            self.spectra[sig['n']][sig['l']][sig['S2p1']] = \
               DummyCXShellSpectrum(self.Z, self.z1,sig['n'],sig['l'], s = sig['S2p1'])

    #    print(Cout, self.spectra)
        emissivity += Cout * self.spectra[sig['n']][sig['l']][sig['S2p1']].calc_spectrum(self.ebins, self.ebins_checksum)


    # return multiplied by velocity
    return emissivity * collvelocity

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

    self.linedataindexes={} # to store the location of all the HDU data


    i = numpy.where((linedata[1].data['Z'] == self.Z) &\
                    (linedata[1].data['z1'] == self.z1))[0]


    nlist = pyatomdb.util.unique(linedata[1].data['n'][i])
    nlist.sort()
    for n in nlist:
       self.linedataindexes[n] = {}

       ii = numpy.where((linedata[1].data['Z'] == self.Z) &\
                        (linedata[1].data['z1'] == self.z1) &\
                        (linedata[1].data['n'] == n))[0]
       if len(ii) > 0:
         llist = pyatomdb.util.unique(linedata[1].data['l'][ii])
         for l in llist:
           self.linedataindexes[n][l] = numpy.where((linedata[1].data['Z'] == self.Z) &\
                        (linedata[1].data['z1'] == self.z1) &\
                        (linedata[1].data['n'] == n) &\
                        (linedata[1].data['l'] == l))[0][0]+2



    self.contdataindexes={} # to store the location of all the HDU data


    i = numpy.where((contdata[1].data['Z'] == self.Z) &\
                    (contdata[1].data['z1'] == self.z1))[0]


    nlist = pyatomdb.util.unique(contdata[1].data['n'][i])
    nlist.sort()
    for n in nlist:
       self.contdataindexes[n] = {}

       ii = numpy.where((contdata[1].data['Z'] == self.Z) &\
                        (contdata[1].data['z1'] == self.z1) &\
                        (contdata[1].data['n'] == n))[0]
       if len(ii) > 0:
         llist = pyatomdb.util.unique(contdata[1].data['l'][ii])
         for l in llist:
           self.contdataindexes[n][l] = numpy.where((contdata[1].data['Z'] == self.Z) &\
                        (contdata[1].data['z1'] == self.z1) &\
                        (contdata[1].data['n'] == n) &\
                        (contdata[1].data['l'] == l))[0][0]+2

    self.linedata = linedata
    self.contdata = contdata

    self.spectra ={}


    self.ebins=0.0
    self.ebins_checksum="0"



  def calc_spectrum(self, collenergy, collvelocity):
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

    emissivity = numpy.zeros(len(self.ebins)-1, dtype=float)
#    ebinshash = hashlib.md5(self.ebins).hexdigest() # tracks for changes to ebins

    # get the l
    l = -1 * (self.acxmodel%4)
    # if l = -4, this returns 0. Correct
    if l==0: l=-4


    for sig in self.crosssectiondata:

      Cout = loginterp(collenergy, sig['E'], sig['C'])

      if Cout > 0:

        if not sig['n'] in self.spectra.keys():
          self.spectra[sig['n']] = {}
        if not l in self.spectra[sig['n']].keys():
          try:
            self.spectra[sig['n']][l] =\
               CXShellSpectrum(self.Z, self.z1, sig['n'], l, \
               self.linedata[self.linedataindexes[sig['n']][l]].data,\
               self.contdata[self.contdataindexes[sig['n']][l]].data)
          except KeyError:
          # Case where there is no data for this Z, z1, nn, l
            self.spectra[sig['n']][l] = \
               DummyCXShellSpectrum(self.Z, self.z1,sig['n'],l)

        emissivity += Cout * self.spectra[sig['n']][l].calc_spectrum(self.ebins, self.ebins_checksum)


    # return multiplied by velocity
    return emissivity * collvelocity

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

    # ok. So. Sort this.
    E_all = numpy.append(Econt_in, eedges)

    cont_tmp = numpy.interp(eedges, Econt_in, cont_in)
    C_all = numpy.append(cont_in, cont_tmp)

    iord = numpy.argsort(E_all)

    # order the arrays
    E_all = E_all[iord]
    C_all = C_all[iord]

    ihi = numpy.where(iord>=n)[0]
    cum_cont = scipy.integrate.cumtrapz(C_all, E_all, initial=0)
    C_out = numpy.zeros(len(eedges))
  #  for i in range(len(eedges)):
  #    arg  = numpy.argwhere(iord==n+i)[0][0]
  #    C_out[i] = cum_cont[ihi[i]]

    C_out = cum_cont[ihi]

    cont = C_out[1:]-C_out[:-1]
    return cont





  def calc_spectrum(self, ebins, ebins_checksum):
    import scipy,scipy.integrate
#    print("WHEE calc_spectrum")

    if self.ebins_checksum == ebins_checksum:
#      print("WHEE no ebins changed")
      return self.spectrum
    else:
      # line emission
#      print("WHEE calculating some stuff!")

      spec = numpy.zeros(len(ebins)-1, dtype=float)

      if len(self.linedata) > 0:
#        print("linedata:", self.linedata)
#        print("contdata:", self.contdata)
#        print("WHEE found %i lines!"%(len(self.linedata)))

        en = pyatomdb.const.HC_IN_KEV_A/self.linedata['Lambda']

        a,b= numpy.histogram(en, bins = ebins, weights = self.linedata['Epsilon'])

        spec += a
      if len(self.contdata) > 0:
        if self.contdata['N_cont'][0] > 2:

          spec += self.expand_E_grid(ebins, self.contdata['E_Cont'][0], self.contdata['Continuum'][0])
#          a = scipy.integrate.cumtrapz(self.contdata['Continuum'][0], x=self.contdata['E_cont'][0], initial=0)
#          b = numpy.interp(ebins,self.contdata['E_cont'][0],a)
#          spec += b[1:]-b[:-1]
        if self.contdata['N_Pseudo'][0] > 2:
          spec += self.expand_E_grid(ebins, self.contdata['E_Pseudo'][0], self.contdata['Pseudo'][0])

#          a = scipy.integrate.cumtrapz(self.contdata['Pseudo'][0], x=self.contdata['E_Pseudo'][0], initial=0)
#          b = numpy.interp(ebins,self.contdata['E_Pseudo'][0],a)
#          spec += b[1:]-b[:-1]
      self.ebins_checksum = ebins_checksum
      self.spectrum = spec
      self.spectrum_ready = True

      return self.spectrum



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


  def calc_spectrum(self, ebins, ebins_checksum):
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

def test():
  # this is a test routine. Let's see what we can do!

  myacx = ACXModel()

#  print('woo1')
#  zzz=input()
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






