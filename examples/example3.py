# List CX line intensities

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

# get the intensity of the O6+ 7->1 (resonance) transition at 1kev/u

collvals = [0.01,0.03,0.1,0.3,1.0,3.0]

w = []
x = []
y = []
z = []

for ic in collvals:
  c = am.calc_line_emissivity(8, 7, 7, 1, collvalue=ic)
  w.append(c['Epsilon'])
  c = am.calc_line_emissivity(8, 7, 6, 1, collvalue=ic)
  x.append(c['Epsilon'])
  c = am.calc_line_emissivity(8, 7, 5, 1, collvalue=ic)
  y.append(c['Epsilon'])
  c = am.calc_line_emissivity(8, 7, 2, 1, collvalue=ic)
  z.append(c['Epsilon'])

w = numpy.array(w)
x = numpy.array(x)
y = numpy.array(y)
z = numpy.array(z)


# make a figure
fig = plt.figure()
fig.show()
ax = fig.add_subplot(211)
ax.loglog(collvals, w, label='resonance')
ax.loglog(collvals, x, label='inter 1')
ax.loglog(collvals, y, label='inter 2')
ax.loglog(collvals, z, label='forbidden')
#ax.plot(wvbins, s6a, label='3.0kev/u', drawstyle='steps')
ax.legend(loc=0)
ax.set_xlabel('Collision param (kev/u)')
ax.set_ylabel('Emissivity (ph cm^3 s^-1 A-1)')


ax2 = fig.add_subplot(212, sharex=ax)
ax2.semilogx(collvals, (x+y+z)/w, label='G ratio')
ax2.legend(loc=0)
ax2.set_xlabel('Collision param (kev/u)')
ax2.set_ylabel('G ration (F+I)/R')

plt.draw()
zzz=input("Press enter to continue")
fig.savefig('example3.pdf') 
fig.savefig('example3.png') 
fig.savefig('example3.svg') 
