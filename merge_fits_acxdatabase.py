import pyatomdb, os, glob, numpy

### Code is designed to copy the fits data from multiple files into one
### representing the data for the cross sections of H and He-like interactions.



fitsdirectory = 'acxdatabase'
donor = 'He'

flist = []
  
flist = glob.glob('%s/*_%s_*fits'%(fitsdirectory, donor))
flist.sort()



filelist = numpy.zeros(len(flist), dtype=numpy.dtype(\
                       {'names':['hdu','Z','z1','resn'],\
                          'formats':[int, int, int, '|S3']}))
hdulist = []

ifile= 0

for f in flist:
  tmpdat = pyatomdb.pyfits.open(f)
  print('tmpdat',tmpdat)
  sigma_cx = numpy.array(tmpdat[1].data)
  
  # change the min and max energies for efficiency
  for i in range(len(sigma_cx)):
     egood= sigma_cx[i]['E'][sigma_cx[i]['C']>0]
     if len(egood) < 2:
       sigma_cx[i]['E_min'] = -1
       sigma_cx[i]['E_max'] = -1
     else:
       sigma_cx[i]['E_min'] = min(egood)
       sigma_cx[i]['E_max'] = max(egood)

  # trim unneeded values.
  sigma_cx = sigma_cx[sigma_cx['E_max']>=0]
  s2p1 = numpy.zeros(len(sigma_cx), dtype=int)
  for i in range(len(s2p1)):
    if numpy.isfinite(sigma_cx[i]['s']):
      s2p1[i] = int(2*sigma_cx[i]['s'])+1
    else:
      s2p1[i] = -1
        
  hdux = pyatomdb.pyfits.BinTableHDU.from_columns(pyatomdb.pyfits.ColDefs(
          [pyatomdb.pyfits.Column(name='n',
             format='1J',
             array=sigma_cx['n']),
           pyatomdb.pyfits.Column(name='l',
             format='1J',
             array=sigma_cx['l']),
           pyatomdb.pyfits.Column(name='S2p1',
             format='1J',
             array=s2p1),
           pyatomdb.pyfits.Column(name='CXType',
             format='1J',
             array=sigma_cx['CXType']),
           pyatomdb.pyfits.Column(name='E_min',
             format='1E',
             unit='keV/amu',
             array=sigma_cx['E_min']),
           pyatomdb.pyfits.Column(name='E_max',
             format='1E',
             unit='keV/amu',
             array=sigma_cx['E_max']),
           pyatomdb.pyfits.Column(name='nC',
             format='1J',
             array=sigma_cx['nC']),
           pyatomdb.pyfits.Column(name='E',
             format='%iE'%(sigma_cx['E'].shape[1]),
             unit='keV/amu',
             array=sigma_cx['E']),
           pyatomdb.pyfits.Column(name='C',
             unit='cm^2',
             format='%iE'%(sigma_cx['C'].shape[1]),
             array=sigma_cx['C'])]))

     
  tmpdat[1].data=hdux.data
  tmpdat[1].header['EXTNAME'] = tmpdat[1].header['System']
  hdulist.append(tmpdat[1].copy())
  print('hdulist',hdulist)
  tmpdat.close()
  #fnameparse
  f2 = f.split('/')
  f3 = f2[1].split('_')
  Z = int(f3[0])
  z1 = int(f3[1])
  filelist[ifile]['Z'] = Z
  filelist[ifile]['z1'] = z1
  resn = 'n'
  if len(numpy.where(hdux.data['l']>=0)[0])>0:
    resn='nl'
  if len(numpy.where(hdux.data['S2p1']>=0)[0])>0:
    resn = 'nls'
  filelist[ifile]['resn'] = resn
  ifile += 1
  
  
hdu1 = pyatomdb.pyfits.BinTableHDU.from_columns(pyatomdb.pyfits.ColDefs(
            [pyatomdb.pyfits.Column(name='Z',
               format='1J',
               array=filelist['Z']),
             pyatomdb.pyfits.Column(name='z1',
               format='1J',
               array=filelist['z1']),
             pyatomdb.pyfits.Column(name='resn',
               format='3A',
               array=filelist['resn'])]))
    

hdu1.header['Donor'] = donor
hdu1.header['EXTNAME'] = 'INDEX'

#print(hdulist[0][0])
hdu1.header['DonMass'] = hdulist[-1].header['DONMASS']



hdu0=pyatomdb.pyfits.PrimaryHDU()
hdulistout  = pyatomdb.pyfits.HDUList([hdu0,hdu1])
for h in hdulist:
  hdulistout.append(h)

hdulistout.writeto('acx2_%s_v1_sigma.fits'%(donor), overwrite=True, checksum=True)
#  tmpCHDUlist = pyatomdb.pyfits.HDUList(CHDUlist)
#  tmpCHDUlist.writeto('fitstmp/%s'%(contfilename), overwrite=True, checksum=True)
