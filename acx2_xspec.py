import numpy, xspec
import acx2 as acx2model

### UPDATED FOR VERSION 2.1.0, 2024-03-03

# CHANGE THESE FILE PATHS TO REFLECT YOUR SYSTEM

# The code will automatically download the ACX2 data files into your
# $ATOMDB directory. If you wish to use non-standard files, 
# specify them here


# tmpdir='/export1/projects/ACX2/'
# Hsigmafile  = tmpdir+'acx2_H_v1_sigma.fits'
# Hlinefile   = tmpdir+'acx2_H_v1_line.fits'
# Hcontfile   = tmpdir+'acx2_H_v1_cont.fits'
# Hesigmafile = tmpdir+'acx2_He_v1_sigma.fits'
# Helinefile  = tmpdir+'acx2_He_v1_line.fits'
# Hecontfile  = tmpdir+'acx2_He_v1_cont.fits'




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
            "vbroad      \"km/s\"    0.0 0.0 0.0 10000.0 10000.0 -0.01",
            "tbroad      \"keV\"     0.0 0.0 0.0 100.0 100.0 -0.01",
            "Redshift      \"\"      0.0 -0.999 -0.999 10.0 10.0 -0.01")

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
             "vbroad      \"km/s\"    0.0 0.0 0.0 10000.0 10000.0 -0.01",
             "tbroad      \"keV\"     0.0 0.0 0.0 100.0 100.0 -0.01",
             "Redshift      \"\"      0.0 -0.999 -0.999 10.0 10.0 -0.01")

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
              "vbroad      \"km/s\"    0.0 0.0 0.0 10000.0 10000.0 -0.01",
              "tbroad      \"keV\"     0.0 0.0 0.0 100.0 100.0 -0.01",
              "Redshift      \"\"      0.0 -0.999 -0.999 10.0 10.0 -0.01")

oneacx2Info = ("$element                14",
               "$ion                    13",
               "collnpar   \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
               "$collntype              1",
               "$acxmodel               8",
               "$recombtype             1",
               "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
               "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
               "vbroad      \"km/s\"    0.0 0.0 0.0 10000.0 10000.0 -0.01",
               "tbroad      \"keV\"     0.0 0.0 0.0 100.0 100.0 -0.01",
               "Redshift      \"\"      0.0 -0.999 -0.999 10.0 10.0 -0.01")


def acx2(engs, params, flux):

  """
  ACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See acx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(acx2, acx2Info, 'add')
    # make a model
    m = xspec.Model('acx2')
  """
  offest = 0
  
  nparam = len(params)
  modeltype=False
  if nparam == len(acx2Info)+1:
    modeltype = 'acx2'
  elif nparam == len(vacx2Info)+1:
    modeltype = 'vacx2'
  elif nparam == len(vvacx2Info)+1:
    modeltype = 'vvacx2'
  elif nparam == len(oneacx2Info)+1:
    modeltype = 'oneacx2'
  else:
    print(modeltype, nparam, len(acx2Info), len(vacx2Info), len(vvacx2Info))

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)

  # vacx model has the 14 main elements
  if modeltype in ['acx2', 'vacx2']:
    elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]
  elif modeltype == 'vvacx2':
    elements = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28]
  elif modeltype=='oneacx2':
    elements = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,28]
    Z = int(params[0])
    ion = int(params[1])

    
  offset = 0
  if modeltype=='oneacx2': offset = 1

  if modeltype == 'acx2':
  # set abundance vector
    abund = numpy.array(params[6+offset])
    acx2_acxmodelobject.set_abund(abund, elements=elements)
  elif modeltype == 'vacx2':
    abund = numpy.array(params[6+offset:20+offset])
    acx2_acxmodelobject.set_abund(abund, elements=elements)
  elif modeltype == 'vvacx2':
    abund = numpy.array(params[6+offset:33+offset])
    acx2_acxmodelobject.set_abund(abund, elements=elements)
  elif modeltype == 'oneacx2':
    abund = numpy.array(params[6+offset])
    acx2_acxmodelobject.set_abund(abund, elements=elements)


  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', elements=elements)

    acx2_acxmodelobject.add_donor('He', elements=elements)
  
#  acx2_acxmodelobject.set_abund(abund)
  
  # get redshift
  redshift = float(params[-2])

  # set energy bins, accounting for redshift
  acx2_acxmodelobject.set_ebins(ebins*(1.0+redshift))

  # check temperature
  if modeltype == 'oneacx2':
    ionfrac = {}
    for i in acx2_acxmodelobject.elements:
      ionfrac[i] = numpy.zeros(i+1, dtype=float)
    ionfrac[Z][ion] = 1.0
#    ionfrac[elements[0]] = numpy.zeros(elements[0]+1, dtype=float)
#   ionfrac[elements[0]][ion]=1.0
    acx2_acxmodelobject.set_ionfrac(ionfrac)


  else:
    acx2_acxmodelobject.set_temperature(params[0])

  # acx fallback model
  acx2_acxmodelobject.set_acxmodel(int(params[3+offset]))

  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4+offset]))

  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5+offset], params[5+offset]])

  # set the collision type (1,2,3 or 4)
  cp = int(params[2+offset])
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)

  

  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(params[1+offset], Tbroaden=params[-3], vbroaden=params[-4])

  # Create a normalization correction factor based on velocity
  if cp==1:
    # convert from center of mass energy to center of mass velocity
    cv = numpy.sqrt(4786031.3*params[1+offset]/25.)
  elif cp == 2:
    # already center of mass velocity
    cv = params[1+offset]
  elif cp == 3:
    # convert from donor velocity (assume donor is H, recombining ion is Carbon-12)
    # c.o.m. vel = (m1v1+m2v2)/(m1+m2)
    cv = 1.0 * params[1+offset]/(1.0+12.0)
  elif cp == 4:
    # convert from recombining ion velocity (assume donor is H, reccombining ion is Carbon-12)
    cv = 12.0 * params[1+offset]/(1.0+12.0)

  # return the flux.
  flux[:] = spec*1e10/cv


def vacx2(engs, params, flux):
  return acx2(engs, params, flux)

def vvacx2(engs, params, flux):
  return acx2(engs, params, flux)
  

def oneacx2(engs, params, flux):
  return acx2(engs, params, flux)



#--def vacx2(engs, params, flux):
#--
#--  """
#--  VACX2 model for data
#--
#--  PARAMETERS
#--  ----------
#--  engs : list[float]
#--    The energy bin edges (from xspec)
#--  params : list[float]
#--    The parameter list. See vacx2Info for definition
#--  flux : list[float]
#--    The array to fill with return values
#--
#--  RETURNS
#--  -------
#--  None
#--    Fills out the flux array with photon cm3 s-1 bin-1, x1e10
#--
#--  USAGE
#--  -----
#--    # load the model into XSPEC
#--    xspec.AllModels.addPyMod(vacx2, vacx2Info, 'add')
#--    # make a model
#--    m = xspec.Model('vacx2')
#--  """
#--
#--  # This is the call that will return everything. So set everything!
#--  ebins = numpy.array(engs)
#--
#--  # vacx model has the 14 main elements
#--  elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]
#--
#--  if len(acx2_acxmodelobject.DonorList)==0:
#--    # Add a donor - hydrogen in this case (default will have H, He)
#--    acx2_acxmodelobject.add_donor('H', \
#--                                  Hlinefile, \
#--                                  Hcontfile, \
#--                                  Hsigmafile,\
#--                                  elements = elements)
#--
#--    acx2_acxmodelobject.add_donor('He', \
#--                                  Helinefile, \
#--                                  Hecontfile, \
#--                                  Hesigmafile,\
#--                                  elements = elements)
#--
#--  # get redshift
#--  redshift = float(params[20])
#--
#--  # set energy bins, accounting for redshift
#--  acx2_acxmodelobject.set_ebins(ebins*(1.0+redshift))
#--
#--  # check temperature
#--  acx2_acxmodelobject.set_temperature(params[0])
#--
#--  acx2_acxmodelobject.set_acxmodel(int(params[3]))
#--
#--  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
#--  acx2_acxmodelobject.set_recombtype(int(params[4]))
#--
#--
#--  # set abundance vector
#--  abund = numpy.array(params[6:20])
#--  acx2_acxmodelobject.set_abund(abund, elements=elements)
#--  # set the H & He fraction from HeFrac
#--  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])
#--
#--  # set the collision type (1,2,3 or 4)
#--  cp = int(params[2])
#--  if cp== 1:
#--    cpunits = 'kev/amu'
#--  elif cp in [2,3,4]:
#--    cpunits = 'km/s'
#--  else:
#--    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))
#--
#--  acx2_acxmodelobject.set_collisiontype(cp, cpunits)
#--
#--
#--  # get the spectrum
#--  spec = acx2_acxmodelobject.calc_spectrum(params[1])
#--
#--
#--  # Create a normalization correction factor based on velocity
#--  if cp==1:
#--    # convert from center of mass energy to center of mass velocity
#--    cv = numpy.sqrt(4786031.3*params[1]/25.)
#--  elif cp == 2:
#--    # already center of mass velocity
#--    cv = params[1]
#--  elif cp == 3:
#--    # convert from donor velocity (assume donor is H, recombining ion is Carbon-12)
#--    # c.o.m. vel = (m1v1+m2v2)/(m1+m2)
#--    cv = 1.0 * params[1]/(1.0+12.0)
#--  elif cp == 4:
#--    # convert from recombining ion velocity (assume donor is H, reccombining ion is Carbon-12)
#--    cv = 12.0 * params[1]/(1.0+12.0)
#--
#--  # return the flux.
#--  flux[:] = spec*1e10/cv
#--
#--
#--
#--
#--
#--
#--
#--
#--
#--def vvacx2(engs, params, flux):
#--
#--  """
#--  VVACX2 model for data
#--
#--  PARAMETERS
#--  ----------
#--  engs : list[float]
#--    The energy bin edges (from xspec)
#--  params : list[float]
#--    The parameter list. See vvacx2Info for definition
#--  flux : list[float]
#--    The array to fill with return values
#--
#--  RETURNS
#--  -------
#--  None
#--    Fills out the flux array with photon cm3 s-1 bin-1, x1e10
#--
#--  USAGE
#--  -----
#--    # load the model into XSPEC
#--    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
#--    # make a model
#--    m = xspec.Model('vvacx2')
#--  """
#--
#--  # This is the call that will return everything. So set everything!
#--  ebins = numpy.array(engs)
#--
#--  # vacx model has the all elements up to nickel except cobalt
#--  elements =[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,\
#--             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,\
#--             21, 22, 23, 24, 25, 26, 28]
#--
#--
#--  if len(acx2_acxmodelobject.DonorList)==0:
#--    # Add a donor - hydrogen in this case (default will have H, He)
#--    acx2_acxmodelobject.add_donor('H', \
#--                                  Hlinefile, \
#--                                  Hcontfile, \
#--                                  Hsigmafile,\
#--                        elements = elements)
#--    acx2_acxmodelobject.add_donor('He', \
#--                                  Helinefile, \
#--                                  Hecontfile, \
#--                                  Hesigmafile,\
#--                        elements = elements)
#--
#--
#--
#--
#--
#--
#--
#--  # get redshift
#--  redshift = float(params[33])
#--
#--  # set energy bins, accounting for redshift
#--  acx2_acxmodelobject.set_ebins(ebins*(1.0+redshift))
#--
#--  # check temperature
#--  acx2_acxmodelobject.set_temperature(params[0])
#--
#--  acx2_acxmodelobject.set_acxmodel(int(params[3]))
#--
#--
#--
#--  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
#--  acx2_acxmodelobject.set_recombtype(int(params[4]))
#--
#--
#--  # set abundance vector
#--  abund = numpy.array(params[6:33])
#--  acx2_acxmodelobject.set_abund(abund, elements=elements)
#--
#--  # set the H & He fraction from HeFrac
#--  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])
#--
#--  # set the collision type (1,2,3 or 4)
#--  cp = int(params[2])
#--  if cp== 1:
#--    cpunits = 'kev/amu'
#--  elif cp in [2,3,4]:
#--    cpunits = 'km/s'
#--  else:
#--    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))
#--
#--  acx2_acxmodelobject.set_collisiontype(cp, cpunits)
#--
#--  # get the spectrum
#--  spec = acx2_acxmodelobject.calc_spectrum(params[1])
#--
#--  # Create a normalization correction factor based on velocity
#--  if cp==1:
#--    # convert from center of mass energy to center of mass velocity
#--    cv = numpy.sqrt(4786031.3*params[1]/25.)
#--  elif cp == 2:
#--    # already center of mass velocity
#--    cv = params[1]
#--  elif cp == 3:
#--    # convert from donor velocity (assume donor is H, recombining ion is Carbon-12)
#--    # c.o.m. vel = (m1v1+m2v2)/(m1+m2)
#--    cv = 1.0 * params[1]/(1.0+12.0)
#--  elif cp == 4:
#--    # convert from recombining ion velocity (assume donor is H, reccombining ion is Carbon-12)
#--    cv = 12.0 * params[1]/(1.0+12.0)
#--
#--  # return the flux.
#--  flux[:] = spec*1e10/cv
#--

# this is how to import the models into pyxspec.
xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
xspec.AllModels.addPyMod(vacx2, vacx2Info, 'add')
xspec.AllModels.addPyMod(acx2, acx2Info, 'add')
xspec.AllModels.addPyMod(oneacx2, oneacx2Info, 'add')

print("Added models acx2, vacx2, vvacx2 and oneacx2 to XSPEC")

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
