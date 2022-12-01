import xspec, numpy

import acx2_xspec

m = xspec.Model('pow+acx2')


# change model parameters however you want here

m.acx2.temperature=0.5
m.acx2.recombtype=2

#####
# bonus feature, new for v1.1.1:
####
# get the emissivity of a line, in this case Ne IX, 7->1, at 1.0keV/amu

l = acx2_xspec.acx2_acxmodelobject.calc_line_emissivity(1.0,10,9,7,1) 

print(l)

# note that this emissivity doesn't take into account any other model effects - 
# e.g. absorption. It's literally just what ACX2 is returning as the
# emissivity of that specific line.

