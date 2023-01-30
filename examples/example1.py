# Create a charge exchange spectrum for a H-He plasma

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
wvbins = numpy.linspace(14.5,19.5,51) # 0.01 A bins
ebins = 12.398425/wvbins[::-1]
am.set_ebins(ebins)


# calculate a spectrum at several different velocities
s1 =am.calc_spectrum(0.01)
s2 =am.calc_spectrum(0.03)
s3 =am.calc_spectrum(0.1)
s4 =am.calc_spectrum(0.3)
s5 =am.calc_spectrum(1.0)
s6 =am.calc_spectrum(3.0)

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
s6a = s6[::-1]/binwidth
s6a = numpy.append(s6a[0], s6a)


# make a figure
fig = plt.figure()
fig.show()
ax = fig.add_subplot(211)
ax.plot(wvbins, s1a, label='0.01kev/u', drawstyle='steps')
ax.plot(wvbins, s2a, label='0.03kev/u', drawstyle='steps')
ax.plot(wvbins, s3a, label='0.1kev/u', drawstyle='steps')
ax.plot(wvbins, s4a, label='0.3kev/u', drawstyle='steps')
ax.plot(wvbins, s5a, label='1.0kev/u', drawstyle='steps')
ax.plot(wvbins, s6a, label='3.0kev/u', drawstyle='steps')
ax.legend(loc=0)
ax.set_xlabel('Wavelength (A)')
ax.set_ylabel('Emissivity (ph cm^3 s^-1 A-1)')


ax2 = fig.add_subplot(212, sharex=ax)
ax2.plot(wvbins, s1a/max(s1a[:len(s1a)//2]), label='0.01kev/u', drawstyle='steps')
ax2.plot(wvbins, s2a/max(s2a[:len(s2a)//2]), label='0.03kev/u', drawstyle='steps')
ax2.plot(wvbins, s3a/max(s3a[:len(s3a)//2]), label='0.1kev/u', drawstyle='steps')
ax2.plot(wvbins, s4a/max(s4a[:len(s4a)//2]), label='0.3kev/u', drawstyle='steps')
ax2.plot(wvbins, s5a/max(s5a[:len(s5a)//2]), label='1.0kev/u', drawstyle='steps')
ax2.plot(wvbins, s6a/max(s6a[:len(s6a)//2]), label='3.0kev/u', drawstyle='steps')
ax2.set_xlabel('Wavelength (A)')
ax2.set_ylabel('Ratio to LymanB')
ax2.legend(loc=0)

 
plt.draw()
zzz=input("Press enter to continue")
fig.savefig('example1.pdf') 
fig.savefig('example1.png') 
fig.savefig('example1.svg') 
