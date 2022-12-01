import acx2, pylab, numpy

# As ever, change these paths as required

tmpdir='/export1/projects/ACX2/'
Hsigmafile  = tmpdir+'acx2_H_v1_sigma.fits'
Hlinefile   = tmpdir+'acx2_H_v1_line.fits'
Hcontfile   = tmpdir+'acx2_H_v1_cont.fits'
Hesigmafile = tmpdir+'acx2_He_v1_sigma.fits'
Helinefile  = tmpdir+'acx2_He_v1_line.fits'
Hecontfile  = tmpdir+'acx2_He_v1_cont.fits'

# For plotting, not necessary outside of this example

fig = pylab.figure(111)
fig.show()
ax = fig.add_subplot(111)

# declare the modelobject
acx2_acxmodelobject = acx2.ACXModel()

# I'm just using oxygen for speed
elements = [8]

# declare wavelength bins
wvbins = numpy.linspace(19,26,71)

# convert to energy bins
ebins = 12.398425/wvbins[::-1]
# set the energy bins                              
acx2_acxmodelobject.set_ebins(ebins)

# add the H and He donors
acx2_acxmodelobject.add_donor('H', \
                              Hlinefile, \
                              Hcontfile, \
                              Hsigmafile,\
                              elements = numpy.array(elements))
acx2_acxmodelobject.add_donor('He', \
                              Helinefile, \
                              Hecontfile, \
                              Hesigmafile,\
                              elements = numpy.array(elements))
                              
# set fallback ACX model and the recombination type (single)
acx2_acxmodelobject.set_acxmodel(8)
acx2_acxmodelobject.set_recombtype(1)

# set relative abundance of H and He donors
acx2_acxmodelobject.set_donorabund(['H','He'], [0.9, 0.1])

# set temperature for ionization balance
acx2_acxmodelobject.set_temperature(0.7)

# set collision parameter to be energy per unit mass
acx2_acxmodelobject.set_collisiontype(1, 'kev/amu')

# return an entire spectrum
s = acx2_acxmodelobject.calc_spectrum(100)

# plot this
ax.plot(ebins, numpy.append(s[0], s), drawstyle='steps')
pylab.draw()

# Get the emissivity of the 4 lines in the He-like O triplet
out = []
for up in [2, 5,6,7]:
  out.append(acx2_acxmodelobject.calc_line_emissivity(100,8,7,up,1))

# print out the emissivities
for o in out:
  print(o, 12.398425/o['Lambda'])
zzz=input()
