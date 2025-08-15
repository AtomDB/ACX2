import acx2, numpy, os

### NOTE CHANGE THIS FILEMAP LOCATION
filemap=os.path.expandvars('$ATOMDB/filemap-acx2_v2.1.0') # the location of the filemap (tells it which files to load)


Z = 26 # Iron
z1 = 17 # The recombining ion charge


specrange = [0.1,0.12] # spectral range in which lines will be listed, in keV
colltype = 1 # energy
collunit = 'kev/amu' # units
collparam= 1.0 # 1kev/u colln

# Iron only ACX2
a = acx2.ACXModel(elements=[26])

# He as only donor
a.add_donor('He')

# Set abundance (unneccesary as there is only 1)
a.set_donorabund(['He'], [1.0])

#set up the ion fraction so it's all 1 ion
z = numpy.zeros(Z+1)
ionfrac = {Z:z}
# This is the pre-recombination balance so [26][17] means
# He + Fe17+ -> He+  + Fe16+
ionfrac[Z][z1]=1.0
a.set_ionfrac(ionfrac)

# set up the collision
a.set_collisiontype(colltype,collunit)

# get linelist
c = a.calc_linelist(collparam, specrange=specrange, specunit='keV')

# sort by emissivity to find strongest lines quickly.
c.sort(order='EPSILON')

# get the full line info (configurations etc)
d= acx2.get_full_line_info(c, filemap = filemap)

print(d) # print out the line list
print(d.dtype) # jsut so you know what's in each column
