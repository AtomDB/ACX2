# Line broadening example

#

import acx2
import matplotlib.pyplot as plt
import numpy

# declare a acx2 model
am = acx2.ACXModel()

# add the H and He donors

am.add_donor('H',
             '/export1/projects/ACX2/acx2_H_v1_line.fits',
             '/export1/projects/ACX2/acx2_H_v1_cont.fits',
             '/export1/projects/ACX2/acx2_H_v1_sigma.fits',
             abundset='AG89',\
             elements=[1,2,6,7,8,10,12,13,14,16,18,20,26,28])
             
am.add_donor('He',
             '/export1/projects/ACX2/acx2_He_v1_line.fits',
             '/export1/projects/ACX2/acx2_He_v1_cont.fits',
             '/export1/projects/ACX2/acx2_He_v1_sigma.fits',
             abundset='AG89',\
             elements=[1,2,6,7,8,10,12,13,14,16,18,20,26,28])


# set the relative abundance of the 2 donors
am.set_donorabund(['H','He'], [0.9, 0.1])

# set the acxmodel to fall back on if there is no other data
am.set_acxmodel(8) # the separable model

# set the collision parameter units
am.set_collisiontype(1, 'kev/u') # keV/amu

# set the recombination type
am.set_recombtype(1) # single recombination. 2= all the way to neutral

# set the temperature
am.set_temperature(1.0) # in keV. determines the ionization balance

# set the binning for the spectrum
wvbins = numpy.linspace(14.5,19.5,1001) # 0.01 A bins
ebins = 12.398425/wvbins[::-1]
am.set_ebins(ebins)

print("This proceses is slow, it will take a while. Sorry")
# calculate a spectrum at several different broadening settings
s1 = am.calc_spectrum(1.0)

print("1/5 spectra done")

# set temperature broadening on 
am.set_broadening(True, broaden_limit=1e-15)
s2 = am.calc_spectrum(1.0)
print("2/5 spectra done")

# set velocity broadening on 
am.set_broadening(True, broaden_limit=1e-15, 
                        velocity_broadening=300.0,\
                        velocity_broadening_units='km/s')
s3 = am.calc_spectrum(1.0)
print("3/5 spectra done")


am.set_broadening(True, broaden_limit=1e-15, 
                        velocity_broadening=600.0,\
                           velocity_broadening_units='km/s')
s4 = am.calc_spectrum(1.0)
print("4/5 spectra done")

am.set_broadening(True, broaden_limit=1e-15, 
                        velocity_broadening=0.0,\
                           thermal_broaden_temperature=5.0)
s5 = am.calc_spectrum(1.0)
print("5/5 spectra done")

# flip the fluxes back into wavelength order, normalize to "per angstrom"
# then prepend a value to it matches the length of ebins (just for plotting)
binwidth = wvbins[1:]-wvbins[:-1]
s1a = s1[::-1]/binwidth
s1a = numpy.append(s1a[0], s1a)
s2a = s2[::-1]/binwidth
s2a = numpy.append(s2a[0], s2a)
s3a = s3[::-1]/binwidth
s3a = numpy.append(s3a[0], s3a)
s4a = s4[::-1]/binwidth
s4a = numpy.append(s4a[0], s4a)
s5a = s5[::-1]/binwidth
s5a = numpy.append(s5a[0], s5a)


# make a figure
fig = plt.figure()
fig.show()
ax = fig.add_subplot(111)
ax.plot(wvbins, s1a, label='none', drawstyle='steps')
ax.plot(wvbins, s2a, label='T=1keV', drawstyle='steps')
ax.plot(wvbins, s3a, label='v=300km/s', drawstyle='steps')
ax.plot(wvbins, s4a, label='v=600km/s', drawstyle='steps')
ax.plot(wvbins, s5a, label='T=5.0keV', drawstyle='steps')
#ax.plot(wvbins, s6a, label='3.0kev/u', drawstyle='steps')
ax.legend(loc=0)
ax.set_xlabel('Wavelength (A)')
ax.set_ylabel('Emissivity (ph cm^3 s^-1 A-1)')

ax.set_xlim(18.85, 19.1)
plt.draw()
zzz=input("Press enter to continue")
fig.savefig('example2.pdf') 
fig.savefig('example2.png') 
fig.savefig('example2.svg') 
