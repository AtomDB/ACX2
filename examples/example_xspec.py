import xspec, numpy

import acx2_xspec

import pylab #(Just for test plots)

m = xspec.Model('pow+acx2')


# change model parameters however you want here

m.acx2.temperature=0.5
m.acx2.recombtype=2

#####
# bonus feature, new for v1.1.1:
####
# get the emissivity of a line, in this case Ne IX, 7->1, at 1.0keV/amu

# to get the emissivity of the current model:
print(acx2_xspec.acx2_acxmodelobject.calc_line_emissivity(10,9,7,1))

# to get the emissivity at other values of the collision parameter:
res={}
vlist=[0.01,0.1,1.0,10.0,100.0, 1000.0]
for up in [2,5,6,7]:
  res[up]=[]
  for vel in vlist:
    m.acx2.collnpar=vel
    l = acx2_xspec.acx2_acxmodelobject.calc_line_emissivity(10, 9, up, 1,\
                                                            collvalue=vel) 
    res[up].append(l['Epsilon'])
    
    print(vel, l)
  
fig = pylab.figure()
fig.show()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax1)
for up in [2,5,6,7]:
  ax1.loglog(vlist, res[up], label=repr(up))
  
ax2.loglog(vlist, (numpy.array(res[7])+1e-30)/(numpy.array(res[2])+1e-30), label='r/f')
ax1.legend()
ax2.legend()
pylab.draw()
zzz=input('Press enter to quit')
# note that this emissivity doesn't take into account any other model effects - 
# e.g. absorption. It's literally just what ACX2 is returning as the
# emissivity of that specific line.

