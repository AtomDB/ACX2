import numpy, xspec
import acx2 as acx2model

# CHANGE THESE FILE PATHS TO REFLECT YOUR SYSTEM
tmpdir='/export1/projects/ACX2/'
Hsigmafile  = tmpdir+'acx2_H_v1_sigma.fits'
Hlinefile   = tmpdir+'acx2_H_v1_line.fits'
Hcontfile   = tmpdir+'acx2_H_v1_cont.fits'
Hesigmafile = tmpdir+'acx2_He_v1_sigma.fits'
Helinefile  = tmpdir+'acx2_He_v1_line.fits'
Hecontfile  = tmpdir+'acx2_He_v1_cont.fits'




#initialize CX object
acx2_acxmodelobject = acx2model.ACXModel()

# These are the definitions XSPEC uses for the inputs

acx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "collnpar   \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
            "$collntype              1",
            "$acxmodel               8",
            "$recombtype             1",
            "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
            "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
            "*redshift     \"\"      0.0")

vacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
             "collnpar    \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
             "$collntype              1",
             "$acxmodel               8",
             "$recombtype             1",
             "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
             "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "*redshift     \"\"      0.0")

vvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
              "collnpar    \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
              "$collntype              1",
              "$acxmodel               8",
              "$recombtype             1",
              "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
              "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Li            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Be            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "B             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "F             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Na            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "P             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cl            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "K             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Sc            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ti            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "V             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cr            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mn            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "*redshift     \"\"      0.0")


bacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "collnpar   \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
            "$collntype              1",
            "$acxmodel               8",
            "$recombtype             1",
            "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
            "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
            "*redshift     \"\"      0.0",
            "$Tbroadening            0",
            "vbroadening   \"km/s\"  0.0 0.01 0.0 0.0 5000.0 5000.")

bvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
             "collnpar    \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
             "$collntype              1",
             "$acxmodel               8",
             "$recombtype             1",
             "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
             "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
             "*redshift     \"\"      0.0",
             "$Tbroadening            0",
             "vbroadening   \"km/s\"  0.0 0.01 0.0 0.0 5000.0 5000.")

bvvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
              "collnpar    \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
              "$collntype              1",
              "$acxmodel               8",
              "$recombtype             1",
              "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
              "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Li            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Be            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "B             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "F             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Na            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "P             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cl            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "K             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Sc            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ti            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "V             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cr            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mn            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "*redshift     \"\"      0.0",
              "$Tbroadening            0",
              "vbroadening   \"km/s\"  0.0 0.01 0.0 0.0 5000.0 5000.")



def acx2(engs, params, flux):

  """
  VACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6]
  redshift = float(params[-2])
  Tbroadening = 0
  vbroadening = 0.0
  abund = numpy.zeros(27)
  for i, j in enumerate([1,2,6,7,8,10,12,13,14,16,18,20,26,27]):
    abund[j-1] = abund_tmp

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec



def bacx2(engs, params, flux):

  """
  BVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6]
  redshift = float(params[-4])
  Tbroadening = int(params[-3])
  vbroadening = float(params[-2])

  abund = numpy.zeros(27)
  for i, j in enumerate([1,2,6,7,8,10,12,13,14,16,18,20,26,27]):
    abund[j-1] = abund_tmp

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec









def vacx2(engs, params, flux):

  """
  VACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6:len(params)-2]
  redshift = float(params[-2])
  Tbroadening = 0
  vbroadening = 0.0
  abund = numpy.zeros(27)
  for i, j in enumerate([1,2,6,7,8,10,12,13,14,16,18,20,26,27]):
    abund[j-1] = abund_tmp[i]

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec



def bvacx2(engs, params, flux):

  """
  BVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6:len(params)-4]
  redshift = float(params[-4])
  Tbroadening = int(params[-3])
  vbroadening = float(params[-2])

  abund = numpy.zeros(27)
  for i, j in enumerate([1,2,6,7,8,10,12,13,14,16,18,20,26,27]):
    abund[j-1] = abund_tmp[i]

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec







def vvacx2(engs, params, flux):

  """
  VVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6:len(params)-2]
  redshift = float(params[-2])
  Tbroadening = 0
  vbroadening = 0.0
  abund = numpy.zeros(27)
  abund[:] = abund_tmp

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec



def bvvacx2(engs, params, flux):

  """
  VVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
  temperature = float(params[0])
  collnpar = float(params[1])
  collntype = int(params[2])
  acxmodel = int(params[3])
  recombtype = int(params[4])
  Hefrac = float(params[5])
  abund_tmp = params[6:len(params)-4]
  redshift = float(params[-4])
  Tbroadening = int(params[-3])
  vbroadening = float(params[-2])
  abund = numpy.zeros(27)
  abund[:] = abund_tmp

  spec = main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening)
  # return the flux.
  flux[:] = spec




def main_acx2(ebins, temperature, collnpar, collntype, acxmodel, recombtype, Hefrac, abund, redshift, Tbroadening, vbroadening):

  """
  VVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!

  # vacx model has the all elements up to nickel except cobalt
  elements =[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,\
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,\
             21, 22, 23, 24, 25, 26, 28]


  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  Hlinefile, \
                                  Hcontfile, \
                                  Hsigmafile,\
                        elements = elements)
    acx2_acxmodelobject.add_donor('He', \
                                  Helinefile, \
                                  Hecontfile, \
                                  Hesigmafile,\
                        elements = elements)

  # set energy bins, accounting for redshift
  acx2_acxmodelobject.set_ebins(ebins*(1.0+redshift))

  # check temperature
  acx2_acxmodelobject.set_temperature(temperature)

  acx2_acxmodelobject.set_acxmodel(int(acxmodel))



  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(recombtype))


  # set abundance vector
  #abund = numpy.array(params[6:33])
  acx2_acxmodelobject.set_abund(abund, elements=elements)

  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-Hefrac, Hefrac])

  # set the collision type (1,2,3 or 4)
  cp = int(collntype)
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)

  acx2_acxmodelobject.set_broadening(Tbroadening, broaden_limit=1e-30,\
                           velocity_broadening=vbroadening, \
                           velocity_broadening_units='km/s')
  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(collnpar)

  # Create a normalization correction factor based on velocity
  if cp==1:
    # convert from center of mass energy to center of mass velocity
    cv = numpy.sqrt(4786031.3*collnpar/25.)
  elif cp == 2:
    # already center of mass velocity
    cv = collnpar
  elif cp == 3:
    # convert from donor velocity (assume donor is H, recombining ion is Carbon-12)
    # c.o.m. vel = (m1v1+m2v2)/(m1+m2)
    cv = 1.0 * collnpar/(1.0+12.0)
  elif cp == 4:
    # convert from recombining ion velocity (assume donor is H, reccombining ion is Carbon-12)
    cv = 12.0 * collnpar/(1.0+12.0)

  # return the flux.
  return(spec*1e10/cv)




# def main_acx2(ebins, redshift, flux):

  # """
  # BVVACX2 model for data

  # PARAMETERS
  # ----------
  # engs : list[float]
    # The energy bin edges (from xspec)
  # params : list[float]
    # The parameter list. See bvvacx2Info for definition
  # flux : list[float]
    # The array to fill with return values

  # RETURNS
  # -------
  # None
    # Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  # USAGE
  # -----
    # # load the model into XSPEC
    # xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # # make a model
    # m = xspec.Model('vvacx2')
  # """

  # # This is the call that will return everything. So set everything!
  # ebins = numpy.array(engs)

  # # vacx model has the all elements up to nickel except cobalt
  # elements =[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,\
             # 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,\
             # 21, 22, 23, 24, 25, 26, 28]


  # if len(acx2_acxmodelobject.DonorList)==0:
    # # Add a donor - hydrogen in this case (default will have H, He)
    # acx2_acxmodelobject.add_donor('H', \
                                  # Hlinefile, \
                                  # Hcontfile, \
                                  # Hsigmafile,\
                        # elements = elements)
    # acx2_acxmodelobject.add_donor('He', \
                                  # Helinefile, \
                                  # Hecontfile, \
                                  # Hesigmafile,\
                        # elements = elements)

  # # set up the abundances
  # if len(params) == 8:
    # # plain acx model
    # abund = numpy.zeros(27)
    # for ind in [1,2,6,7,8,10,12,13,14,16,18,20,26,27]:
      # abund[ind-1] = params[6]
    # DoTBroadening = False
    # v_broadening = 0.0

  # elif len(params) == 10:
    # # bacx model
    # abund = numpy.zeros(27)
    # for ind in [1,2,6,7,8,10,12,13,14,16,18,20,26,27]:
      # abund[ind-1] = params[6]
    # DoTBroadening = params[-2]
    # v_broadening = params[-1]

  # elif len(params) == 21:
    # # vacx model
    # abund = numpy.zeros(27)
    # for i, ind in enumerate([1,2,6,7,8,10,12,13,14,16,18,20,26,27]):
      # abund[ind-1] = params[6+i]
    # DoTBroadening = False
    # v_broadening = 0.0
  # elif len(params) == 23:
    # # bvacx model
    # abund = numpy.zeros(29)
    # for ind in [1,2,6,7,8,10,12,13,14,16,18,20,26,27]:
      # abund[ind-1] = params[6]
    # DoTBroadening = params[-2]
# #    v_broadening = params[-1]
# #
# #  elif len(params) == 34:
# #    # vvacx model
# #    abund = numpy.zeros(27)
# #    for i, ind in enumerate([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
# #                             18,19,20,21,22,23,24,25,26,27]):
# #      abund[ind-1] = params[6+i]
# #    DoTBroadening = False
# #    v_broadening = 0.0
# #  elif len(params) == 36:
# #    # bvacx model
# #    abund = numpy.zeros(29)
# #    for i, ind in enumerate([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
# #                             18,19,20,21,22,23,24,25,26,27]):
# #      abund[ind-1] = params[6+i]
# #    DoTBroadening = params[-2]
# #    v_broadening = params[-1]


  # # get redshift
  # redshift = float(params[-3])

  # # set energy bins, accounting for redshift
  # acx2_acxmodelobject.set_ebins(ebins*(1.0+redshift))

  # # check temperature
  # acx2_acxmodelobject.set_temperature(params[0])

  # acx2_acxmodelobject.set_acxmodel(int(params[3]))



  # # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  # acx2_acxmodelobject.set_recombtype(int(params[4]))


  # # set abundance vector
  # abund = numpy.array(params[6:33])
  # acx2_acxmodelobject.set_abund(abund, elements=elements)

  # # set the H & He fraction from HeFrac
  # acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])

  # # set the collision type (1,2,3 or 4)
  # cp = int(params[2])
  # if cp== 1:
    # cpunits = 'kev/amu'
  # elif cp in [2,3,4]:
    # cpunits = 'km/s'
  # else:
    # raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  # acx2_acxmodelobject.set_collisiontype(cp, cpunits)

  # DoTBroadening = params[34]
  # if DoTBroadening == 0:
    # DoTBroadening = False
  # else:
    # DoTBroadening = True
  # v_broadening = float(params[35])

  # acx2_acxmodelobject.set_broadening(DoTBroadening, broaden_limit=1e-30,\
                           # velocity_broadening=v_broadening, \
                           # velocity_broadening_units='km/s')


  # # get the spectrum
  # spec = acx2_acxmodelobject.calc_spectrum(params[1])

  # # Create a normalization correction factor based on velocity
  # if cp==1:
    # # convert from center of mass energy to center of mass velocity
    # cv = numpy.sqrt(4786031.3*params[1]/25.)
  # elif cp == 2:
    # # already center of mass velocity
    # cv = params[1]
  # elif cp == 3:
    # # convert from donor velocity (assume donor is H, recombining ion is Carbon-12)
    # # c.o.m. vel = (m1v1+m2v2)/(m1+m2)
    # cv = 1.0 * params[1]/(1.0+12.0)
  # elif cp == 4:
    # # convert from recombining ion velocity (assume donor is H, reccombining ion is Carbon-12)
    # cv = 12.0 * params[1]/(1.0+12.0)

  # # return the flux.
  # flux[:] = spec*1e10/cv

# this is how to import the models into pyxspec.
#xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
#xspec.AllModels.addPyMod(vacx2, vacx2Info, 'add')
#xspec.AllModels.addPyMod(acx2, acx2Info, 'add')

#xspec.AllModels.addPyMod(bvvacx2, bvvacx2Info, 'add')
#xspec.AllModels.addPyMod(bvacx2, bvacx2Info, 'add')
#xspec.AllModels.addPyMod(bacx2, bacx2Info, 'add')
xspec.AllModels.addPyMod(NEWvvacx2, bvvacx2Info, 'add')

print("Models acx2, vacx2, vvacx2, bacx2, bvacx2, bvvacx2 succesfully imported")
print("To use, load in the usual way, e.g. xspec.Model('vacx2+pow')")

# m = xspec.Model('vvacx2')
# m = xspec.Model('vacx2')
# m = xspec.Model('acx2')





#fs1 = xspec.FakeitSettings(response='acisi_aimpt_cy20.rmf', \
#                           arf='acisi_aimpt_cy20.arf',\
#                           exposure='1e6', \
#                           fileName='fake_acx2.fak')

#xspec.AllData.fakeit(nSpectra=1, settings=fs1, applyStats=True)

#xspec.Plot.device='/xs'
#xspec.Plot.xAxis='A'

#xspec.Plot('data')
