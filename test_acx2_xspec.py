import acx2_xspec, numpy, pylab


### test 1. Make a spectrum of $O^{7+} + H \rightarrow O^{6+} + H^{+}$


Z=9
m1 = acx2_xspec.acx2model.ACXDonorModel('H', \
                '/export1/projects/ACX2/acx2_H_v1_line.fits', \
                '/export1/projects/ACX2/acx2_H_v1_cont.fits', \
                '/export1/projects/ACX2/acx2_H_v1_sigma.fits', \
                 elements=[6])
#                 elements=[1,2,6,7,8,9,10,12,13,14,16,18,20,26,28])


wvbins = numpy.linspace(10,400, 10001)
ebins = 12.398425/(wvbins[::-1])
m1.set_ebins(ebins)




#m1.set_collisiontype(1, 'kev/amu')

m1.set_collisiontype(2, 'cm/s')
#m1.set_collisiontype(3, 'cm/s')
#m1.set_collisiontype(4, 'cm/s')


m1.set_temperature(0.6)
#m1.ionfrac = {}
#m1.ionfrac[Z] = numpy.zeros(Z+1)

#m1.ionfrac[Z][-2]=1.0

v = 600e5
cp3= v*17
cp4 = v*17/16.

s = m1.calc_spectrum(v) # param2


fig = pylab.figure()
fig.show()
ax = fig.add_subplot(111)

ax.semilogy(12.398425/m1.ebins[::-1], numpy.append(0,s[::-1])+1e-40, drawstyle='steps')

# now change the ion
#m1.ionfrac[Z] = numpy.zeros(Z+1)
#m1.ionfrac[Z][-1]=1.0
#s = m1.calc_spectrum(v) # param2
#ax.semilogy(12.398425/m1.ebins[::-1], numpy.append(0,s[::-1])+1e-40, drawstyle='steps')


#m1.ionfrac[Z] = numpy.zeros(Z+1)
#m1.ionfrac[Z][-3]=1.0
#s = m1.calc_spectrum(v) # param2
#ax.semilogy(12.398425/m1.ebins[::-1], numpy.append(0,s[::-1])+1e-40, drawstyle='steps')

pylab.draw()
zzz=input()
