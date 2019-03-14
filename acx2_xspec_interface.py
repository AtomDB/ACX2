import acx2_xspec, numpy, xspec

#initialize CX object
m1 = acx2_xspec.ACXModel()

# Add a donor - hydrogen in this case (default will have H, He)
m1.add_donor('H', \
                'tmp_cxline2.fits', \
                'tmp_cxcont.fits', \
                'tmp_cx_sigma.fits',\
                elements=[1,2,6,7,8,10,12,13,14,16,18,20,26,28])




acx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862  86.  86.  0.01",
               "collenergy    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01",
               "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
               "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
               "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01")

vacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862  86.  86.  0.01",
               "collenergy    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01",
               "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
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


vvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862  86.  86.  0.01",
               "collenergy    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01",
               "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
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


def vacx(engs, params, flux):

  """
  This is the module to import in to xspec, using a command like:
  xspec.AllModels.addPyMod(vvacx2, vvacxInfo, 'add')

  xspec.Model('vvacx2')
  """
#  input('fuckit1')
  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)
#  input('fuckit2')

  # check energies
  m1.set_ebins(ebins)
#  input('fuckit3')

  # check temperature
  m1.set_temperature(params[0])
#  input('fuckit4')

  m1.set_acxmodel(int(params[2]))
#  input('fuckit5')

  elements =[1,2,6,7,8,10,12,13,14,16,18,20,26,28]
  abund = numpy.array(params[4:18])
  m1.set_abund(elements, abund)
#  m1.set_donorabund(['H','He'], [1-params[3], params[3]])

  m1.set_collisiontype(1, 'kev/amu')

  spec = m1.calc_spectrum(params[1])

  for i in range(len(flux)):
    flux[i]=float(spec[i])*1e10
#    print("flux[%i]=%e"%(i,flux[i]))




def test11(engs, params, flux):
  for i in range(len(engs)-1):
    flux[i] = engs[i]*1


test11Info = ("temperature   \"keV\"   1.0 0.00862 0.00862  86.  86.  0.01",
               "collenergy    \"kev/u\" 1.0 0.01 0.2 100. 1000. 0.01")



xspec.AllModels.addPyMod(vacx, vacx2Info, 'add')
m = xspec.Model('vacx')

#xspec.AllModels.addPyMod(test11, test11Info, 'add')
#m = xspec.Model('test11')

fs1 = xspec.FakeitSettings(response='acisi_aimpt_cy20.rmf', \
                           arf='acisi_aimpt_cy20.arf',\
                           exposure='1e6', \
                           fileName='fake_acx2.fak')

xspec.AllData.fakeit(nSpectra=1, settings=fs1, applyStats=True)

xspec.Plot.device='/xs'
xspec.Plot.xAxis='A'

xspec.Plot('data')
