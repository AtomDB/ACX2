import xspec, numpy

import acx2_xspec

import pylab #(Just for test plots)

m = xspec.Model('pow+vacxnei')


### NOTE: this will use the most recent ionization balance data
### from AtomDB. This is currently v3.1.X
### If this is not the data set you want to use you will need to 
### update your pyatomdb.util.switch_version('3.0.9') (or whatever) to 
### choose an older version.
### A smarter way of doing this is on the to-do list.

# change model parameters however you want here

m.vacxnei.temperature=0.5
m.vacxnei.recombtype=2

# set a recombining plasma
m.vacxnei.kT_init=2.0
m.vacxnei.Tau=1e11

# let's check the ionization fraction
print(acx2_xspec.acx2_acxmodelobject.ionfrac)

# can compare to xspec nei model values if you set chatter=25 in xspec
# xspec.Xset.chatter=25

# then make spectra in the usual way


xspec.Plot.device='/null'
xspec.Plot('model')

x = xspec.Plot.x()
y = xspec.Plot.model()



fig = pylab.figure()
fig.show()
ax1 = fig.add_subplot(111)
ax1.plot(x,y, label='Tau=1e11')


# change the ionization timescale
m.vacxnei.Tau=1e13
xspec.Plot('model')

x = xspec.Plot.x()
y = xspec.Plot.model()
ax1.plot(x,y, label='Tau=1e13')
ax1.legend(loc=0)

ax1.set_xlabel(xspec.Plot.labels()[0])
ax1.set_ylabel(xspec.Plot.labels()[1])
ax1.set_title(xspec.Plot.labels()[2])

  
pylab.draw()
zzz=input('Press enter to quit')
# note that this emissivity doesn't take into account any other model effects - 
# e.g. absorption. It's literally just what ACX2 is returning as the
# emissivity of that specific line.
