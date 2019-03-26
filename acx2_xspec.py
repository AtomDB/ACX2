import acx2_xspec, numpy, xspec

#initialize CX object
acx2_acxmodelobject = acx2_xspec.ACXModel()





# These are the definitions XSPEC uses for the inputs

acx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "collnpar   \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
            "collntype     \"\"      1.0 1.0 1.0 4.0 4.0 -0.01",
            "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
            "recombtype    \"\"      1.0 1.0 1.0 2.0 2.0 -0.01",
            "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
            "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01")

vacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
             "collnpar    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01",
             "collntype     \"\"      1.0 1.0 1.0 4.0 4.0 -0.01",
             "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
             "recombtype    \"\"      1.0 1.0 1.0 2.0 2.0 -0.01",
             "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
             "H             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "He            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "C             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "N             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "O             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "F             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Al            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Si            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "S             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
             "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 0.01")


vvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
              "collnpar    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01",
              "collntype     \"\"      1.0 1.0 1.0 4.0 4.0 -0.01",
              "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
              "recombtype    \"\"      1.0 1.0 1.0 2.0 2.0 -0.01",
              "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
              "H             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "He            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Li            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Be            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "B             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "C             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "N             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "O             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "F             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Na            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Al            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Si            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "P             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "S             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Cl            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "K             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Sc            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Ti            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "V             \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Cr            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Mn            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 0.01",
              "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 0.01")

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

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)

  # vacx model has the 14 main elements
  elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]

  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  'acx2_H_v1_line.fits', \
                                  'acx2_H_v1_cont.fits', \
                                  'acx2_H_v1_sigma.fits',\
                                  elements = numpy.array(elements))

    acx2_acxmodelobject.add_donor('He', \
                                  'acx2_He_v1_line.fits', \
                                  'acx2_He_v1_cont.fits', \
                                  'acx2_He_v1_sigma.fits',\
                                  elements = numpy.array(elements))

  # check energies
  acx2_acxmodelobject.set_ebins(ebins)

  # check temperature
  acx2_acxmodelobject.set_temperature(params[0])

  acx2_acxmodelobject.set_acxmodel(int(params[3]))


  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4]))


  # set abundance vector
  abund = numpy.array(params[6])
  acx2_acxmodelobject.set_abund(abund, elements=elements)

  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])

  # set the collision type (1,2,3 or 4)
  cp = int(params[2])
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)


  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(params[1])

  # return the flux.
  flux[:] = spec*1e10


def vacx2(engs, params, flux):

  """
  VACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vacx2, vacx2Info, 'add')
    # make a model
    m = xspec.Model('vacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)

  # vacx model has the 14 main elements
  elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]

  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  'acx2_H_v1_line.fits', \
                                  'acx2_H_v1_cont.fits', \
                                  'acx2_H_v1_sigma.fits',\
                                  elements = elements)

    acx2_acxmodelobject.add_donor('He', \
                                  'acx2_He_v1_line.fits', \
                                  'acx2_He_v1_cont.fits', \
                                  'acx2_He_v1_sigma.fits',\
                                  elements = elements)

  # check energies
  acx2_acxmodelobject.set_ebins(ebins)

  # check temperature
  acx2_acxmodelobject.set_temperature(params[0])

  acx2_acxmodelobject.set_acxmodel(int(params[3]))

  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4]))


  # set abundance vector
  abund = numpy.array(params[6:20])
  acx2_acxmodelobject.set_abund(abund, elements=elements)
  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])

  # set the collision type (1,2,3 or 4)
  cp = int(params[2])
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)


  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(params[1])

  # return the flux.
  flux[:] = spec*1e10









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

  # vacx model has the all elements up to nickel except cobalt
  elements =[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,\
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,\
             21, 22, 23, 24, 25, 26, 28]


  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  'acx2_H_v1_line.fits', \
                                  'acx2_H_v1_cont.fits', \
                                  'acx2_H_v1_sigma.fits',\
                        elements = elements)
    acx2_acxmodelobject.add_donor('He', \
                                  'acx2_He_v1_line.fits', \
                                  'acx2_He_v1_cont.fits', \
                                  'acx2_He_v1_sigma.fits',\
                        elements = elements)







  # check energies
  acx2_acxmodelobject.set_ebins(ebins)

  # check temperature
  acx2_acxmodelobject.set_temperature(params[0])

  acx2_acxmodelobject.set_acxmodel(int(params[3]))



  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4]))


  # set abundance vector
  abund = numpy.array(params[6:33])
  acx2_acxmodelobject.set_abund(abund, elements=elements)

  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])

  # set the collision type (1,2,3 or 4)
  cp = int(params[2])
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)

  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(params[1])

  # return the flux.
  flux[:] = spec*1e10



# this is how to import the models into pyxspec.
xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
xspec.AllModels.addPyMod(vacx2, vacx2Info, 'add')
xspec.AllModels.addPyMod(acx2, acx2Info, 'add')



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
