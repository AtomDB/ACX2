import pyatomdb, sys, numpy, pickle, time, os, re
###
###
# This files contains all the scripts required to create the various
# CX fits from scratch. This is a big ask! FUN will ensue.
# writte in python 3, for now. 3.5 compatible definitely
###
###

UNIVERSAL_CX_CROSSSECTION    = 1e-10

def calc_l_dstn(acx_ldist, n, z):
  """
  Calculate the l distribution for shell n, l distribution

  PARAMETERS
  ----------
  acx_ldist: string
     "even" : even statistical distribution
     "statistical" : statistical weighting
     "landauzener" : Landau-Zener from Janev & Winter
     "separable" : Separable from Janev & Winter

  n : int
    principal quantum number

  z : int
    recombining ion charge

  RETURNS
  -------
  ldist : array
    Fraction of capture into each l, normalized to 1
  """

  ldist = numpy.zeros(n, dtype=float)
  print("n=",n)
  if acx_ldist == 'even':
    ldist[:]=1.0/n

  elif acx_ldist == 'statistical':
    ldist[:]=numpy.arange(n)*4.0+2
    ldist/=sum(ldist)

  elif acx_ldist == 'landauzener':
    llist = numpy.arange(n, dtype=float)
    ldist = numpy.zeros(n, dtype=float)
    for il,l in enumerate(llist):
      if n-l-1 < 0:
        ldist[il] = 0
      elif n-2 < 0:
        ldist[il] = 0
      else:
        ldist[il]= (l *(l+1) * (2*l+1) * numpy.math.factorial(n-1) * \
               numpy.math.factorial(n-2))*1.0/\
           (numpy.math.factorial(n+l)*numpy.math.factorial(n-l-1))
    if sum(ldist) > 0:
      ldist /=sum(ldist)

  elif acx_ldist == 'separable':
    llist = numpy.arange(n, dtype=int)
    ldist = ( (2*llist+1)*1.0/z) * numpy.exp( (-llist*(llist+1)*1.0/z))
    ldist /=sum(ldist)

  return ldist

def get_ion_sigma_acx(Z,z1, acxtype, donor, datacache=False, ndist = False, nmax=False, nmin=False, debug = False):
  """
  Make the ACX spectrum for a particular ion (no Kronos data)

  PARAMETERS
  ----------

  Z: int
    Atomic No of element
  z1: int
    Recombining ion charge+1
  axctype: int
    From 1-4, for different types:\n
    1: even\n
    2: statistical\n
    3: landauzener\n
    4: separable\n
    See AtomDB CX manual for details
  donor: string
    element symbol for donor ion. Currently "H" or "He" accepted.
  datacache: dict
    pass to the code to speed calculations (caches results)
  ndist : bool
    If True, calculate spectra separately for each n shell. Otherwise
    calculate entire spectrum for a given n'.
  nmax : int
    If ndist is True, this is the maximum n to include
  nmin : int
    If ndist is True, this is the minimum n to include
  debug : bool
    If true, create debugging file


  RETURNS
  -------
  sigma: cross sections, normalized to 3e-15, for each n shell
  """

  settings={}
  settings['filemap'] = 'filemap_v3.0.9_cx_2'
  if datacache==False:
    datacache={}
  if donor == 'H':
    donor_ionpot = pyatomdb.atomdb.get_ionpot(1,1) # in keV
  else:
    donor_ionpot = pyatomdb.atomdb.get_ionpot(2,1) # in keV


  #
  # GREAT!
  # now code up the things!

  # Find the ACX model of interest.

  q = float(z1)-1 # recombining ion charge

  # get recombining ion potential
  ip = pyatomdb.atomdb.get_ionpot(Z,z1-1)

  # set up array of nshells to calculate
  if ndist:
    nlist = numpy.arange(nmin, nmax+1, dtype=int)
    nfrac = numpy.ones(len(nlist), dtype=float)
  else:
    ncap = q * numpy.sqrt(pyatomdb.const.RYDBERG*1000/donor_ionpot) /((1+ (q-1)/numpy.sqrt(2*q))**0.5)
    nlist = numpy.array([int(numpy.floor(ncap)), int(numpy.ceil(ncap))], dtype=int)
    nfrac = numpy.array([1-(ncap-nlist[0])/(nlist[1]-nlist[0]), (ncap-nlist[0])/(nlist[1]-nlist[0]) ], dtype=float)

  shlmax =1
  n = 1
  l=0
  while n < max(nlist)+1:
    l+=1
    if l==n:
      l=0
      n+=1
    shlmax+=1

  #nlo = numpy.floor(ncap)
  #nup =nlo+1
  #energies = numpy.linspace(1,1e6,11)
  #fnup = (ncap-nlo)
  #fnlo = 1-fnup
  #nfrac=[fnlo,fnup]
  #nlev = nlo+nup
  # recombined data
  lvdat = pyatomdb.atomdb.get_data(Z,z1-1,'LV', settings=settings, datacache=datacache)
  ladat = pyatomdb.atomdb.get_data(Z,z1-1,'LA', settings=settings, datacache=datacache)
  if lvdat==False: return numpy.zeros(10)
  # Define each level by it's state, relevant to the distribution here
  lvstate = numpy.zeros(len(lvdat[1].data), dtype=int)
  lvlstate = numpy.zeros(len(lvdat[1].data), dtype=int)
  lvLstate = numpy.zeros(len(lvdat[1].data), dtype=int)

  # First, autoionizing states are RIGHT OUT
  lvstate[lvdat[1].data['AAUT_TOT']>0] = -1 #no good

  # number of electrons in recombined ion
  nel = Z+2-z1
  for i in range(len(lvdat[1].data)):
    tmp,x = pyatomdb.atomic.config_to_occup(lvdat[1].data['ELEC_CONFIG'][i], nel=nel)
    lvdat[1].data['ELEC_CONFIG'][i]=pyatomdb.atomic.occup_to_cfg(tmp)

  # now, look at the parent ion ground state, and find places to add an electron.
  nel_p = nel-1
  if z1 <= Z:
    plvdat = pyatomdb.atomdb.get_data(Z, z1, 'LV',settings=settings, datacache=datacache)
    print("cfg=",plvdat[1].data['ELEC_CONFIG'][0], 'nel=', nel_p, 'shlmax=', shlmax)
    try:
      pgndocc,x = pyatomdb.atomic.config_to_occup(plvdat[1].data['ELEC_CONFIG'][0], nel = nel_p, shlmax = shlmax)
    except IndexError:
      pgndocc,x = pyatomdb.atomic.config_to_occup(plvdat[1].data['ELEC_CONFIG'][0], nel = nel_p)
    print("pgndocc=",pgndocc)

  else:
    pgndocc = numpy.zeros(shlmax, dtype=int)

  n=1
  l=0



  for i in range(len(pgndocc)):
    occup = (l*4)+2
    tmpocc = pgndocc * 1
    if tmpocc[i] < occup:
      # add the electron
      tmpocc[i] += 1

      tmpcfg = pyatomdb.atomic.occup_to_cfg(tmpocc)

      j = numpy.where(lvdat[1].data['ELEC_CONFIG']==tmpcfg)[0]

      if len(j) == 0: continue

      lvlstate[j] = l
      for jj in j:
        if lvstate[jj] != -1:
          lvstate[jj] = n

      lvLstate[j] = lvdat[1].data['L_QUAN'][j]
#        lvstate[j] = lvdat[1].data['N_QUAN'][j]
    l+=1
    if l==n:
      l=0
      n+=1

  acx_ltype='l'
  if acxtype==1:
    acx_ldist='even'
  elif acxtype==2:
    acx_ldist='statistical'
  elif acxtype==3:
    acx_ldist='landauzener'
  elif acxtype==4:
    acx_ldist='separable'

  # Now start putting things in to the correct shells

  TOTSIGMA=3e-15 #cm^-2

  sigma = numpy.zeros(len(lvdat[1].data), dtype=float)

  # find the fraction in each n

  # for each n, find how many levels have a legit capture chance


  nmap ={}
  for ii, n in enumerate(nlist):
#    print(lvstate.shape)
#   print(lvlstate.shape)
    i = numpy.where( (lvstate==n) & \
                     (lvlstate>-1))[0]
    iii=0
    while len(i) == 0:

      print("n = %i has no potential capture levels. Promoting it"%(n+iii))
      iii+=1
      i = numpy.where( (lvstate==n+iii) & \
                       (lvlstate>-1))[0]

    nmap[n] = n+iii


  # find the fraction in each l
  sigma ={}
  # go through each n
  print('nlist:', nlist)
  print('nfrac:', nfrac)

  print('nmap:', nmap)


  for ii, n in enumerate(nlist):
    if n==0: continue
    sigma[n] = numpy.zeros(len(lvstate), dtype=float)
    # calculate the fraction in each L or l
    ntmp = nmap[n]
    ldist = calc_l_dstn(acx_ldist, ntmp, q)
    # find all the levels with the correct n
    nlist = numpy.where(lvstate==ntmp)[0]
    if len(nlist) == 0:
      print("no levels with n=%i here"%(n))
    else:
      # find the range of possible l or L for these levels
#      if acx_ltype=='L':
#        tmplstate = lvLstate*1
#        llist = pyatomdb.util.unique(lvLstate[nlist])
#        # if there are -1s, assign them lvLstate = lvlstate
#        if -1 in llist:
#          tmplstate[lvLstate==-1] = lvlstate[lvLstate==-1]
#        llist = pyatomdb.util.unique(tmplstate[nlist])
      if acx_ltype=='l':
        tmplstate = lvlstate*1

        llist = pyatomdb.util.unique(lvlstate[nlist])

      # find any unmatched l, renormalise it out
      for l in range(ntmp):
        if not l in llist:
          ldist[l]=0
          ldist/=sum(ldist)
      print("sum ldist", sum(ldist), ldist)
      # for each l or L (using tmplstate), split up the sigma
      for l in range(ntmp):
        # go through each n, find possible l
        ilmatch = numpy.where((lvstate==ntmp) & (tmplstate==l))[0]
        if len(ilmatch)==0:
          continue

        # weight by statistical weight
        sigma[n][ilmatch] = nfrac[iii] * TOTSIGMA * ldist[l] * \
                            lvdat[1].data['LEV_DEG'][ilmatch]*1.0/sum(lvdat[1].data['LEV_DEG'][ilmatch])

  if debug:
    f=open('sigma_mode_%i_donor_%s_ion_%i_%i.dat'%(acxtype, donor,Z, z1), 'w')
    s = 'level'
    nn = list(sigma.keys())
    nn.sort()
    for nnn in nn:
      s  = s+ "     n=%2i    "%(nnn)
    s+="\n"

    f.write(s)
    s=repr(nfrac)+'\n'
    f.write(s)
    for i in range(len(lvstate)):
      s = '%5i '%(i+1)
      for nnn in nn:
        s = s+' %e'%(sigma[nnn][i])
      s+=" %s"%(lvdat[1].data['ELEC_CONFIG'][i])
      s+="\n"
      f.write(s)

    f.close()


  return sigma


def create_spectrum(Z, z1_recombining,acxtype, donor, sigma, datacache=False):
  z1=z1_recombining-1
  settings={}
  VELOCITY= 100.0 #1m/s

  settings['filemap'] = 'filemap_v3.0.9_cx_2'
  # some continuum settings
  settings['LinearGrid'] = True
  settings['GridMinimum'] = 0.001
  settings['GridMaximum'] = 100.
  settings['NumGrid'] = 100000
  settings['MinEpsilon']=1e-10
  settings['TwoPhoton']=True


  if datacache==False:
    datacache={}
  te = 1.0
  dens = 1.0

  ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                            settings['GridMinimum'], \
                            settings['GridMaximum'], \
                            settings['NumGrid'])


  matrixA_in = {}
  matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
       pyatomdb.apec.gather_rates(Z, z1, te, dens, datacache=datacache,\
                 settings=settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

  nlev = len(sigma)
  matrixA = numpy.zeros([nlev,nlev],dtype=float)
  continuum = {}

  for i in range(len(matrixA_in['init'])):
    matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]

  if sum(sigma[1:])>0:
    matrixB = -1*sigma*VELOCITY
    print("sum matrixB", sum(matrixB[1:]))
    print("matrix B > 0:",sum(matrixB[1:]>0))
    levpop_this = pyatomdb.apec.calc_cascade_population(matrixA, matrixB)
  else:
    levpop_this = numpy.zeros(nlev)

  if sum(levpop_this) > 0:
    # correct for abundance
    levpop_this *= pyatomdb.atomdb.get_abundance()[Z]
    # get lines!
    linelist, twophot = \
      pyatomdb.apec.do_lines(Z, z1, levpop_this , 1.0, datacache=datacache, settings=settings, z1_drv_in=z1_recombining)
    continuum['twophot']= twophot

    print("Start filtering linelist at %s"%(time.asctime()))

    MinEpsilon = settings['MinEpsilon']
    ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                              settings['GridMinimum'], \
                              settings['GridMaximum'], \
                              settings['NumGrid'])

  #  istrong = linelist>=MinEpsilon
    pseudocont = numpy.zeros(len(ebins)-1, dtype=float)

    if len(linelist) > 0:
      linelist = linelist[(linelist['lambda']>pyatomdb.const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                          (linelist['lambda']<pyatomdb.const.HC_IN_KEV_A /settings['GridMinimum'])]

    if len(linelist) > 0:
      weaklines = linelist[(linelist['epsilon']< MinEpsilon)]

      print("filtered out % i weak lines, left %i strong lines"%(len(weaklines), len(linelist)))

      pseudocont,zzz = numpy.histogram(pyatomdb.const.HC_IN_KEV_A/weaklines['lambda'], bins=ebins, weights=weaklines['epsilon'])
#      for iline, line in enumerate(weaklines):
#        e = pyatomdb.const.HC_IN_KEV_A /line['lambda']
#        ibin = numpy.where(ebins>e)[0][0] - 1
#        if (iline+1)%100000==0:
#          print "  on line %i"%(iline)
#        pseudocont[ibin]+=line['epsilon']

      linelist = linelist[linelist['epsilon'] >= MinEpsilon]
    print("Finish filtering linelist at %s"%(time.asctime()))


    econt, contin = pyatomdb.apec.compress_continuum(ebins, continuum['twophot'], pyatomdb.const.TOLERANCE, minval=1e-38)

    epseudo, pseudo = pyatomdb.apec.compress_continuum(ebins, pseudocont, pyatomdb.const.TOLERANCE, minval=1e-38)

    ret={}
    ret['lines']=linelist
    cont = {}
    cont['E_Cont'] = econt
    cont['Cont'] = contin
    cont['E_Pseudo'] = epseudo
    cont['Pseudo'] = pseudo
    ret['cont'] = cont
    ret['Z'] = Z
    ret['rmJ'] = z1_recombining
    ret['settings'] = settings
    ret['acxtype'] = acxtype
    ret['donor'] = donor
    ret['velocity'] = VELOCITY


    pickle.dump(ret, open('results/acx_%i_donor_%s_Z_%i_rmJ_%i.pkl'%(acxtype,donor,Z,z1_recombining), 'wb'))


def process_ion_acx(Z,z1, acxtype, donor, datacache=False, ndist = False, nmax=False, nmin=False, debug = False):
  """
  Make the ACX spectrum for a particular ion (no Kronos data)

  PARAMETERS
  ----------

  Z: int
    Atomic No of element
  z1: int
    Recombining ion charge+1
  axctype: int
    From 1-4, for different types:\n
    1: even\n
    2: statistical\n
    3: landauzener\n
    4: separable\n
    See AtomDB CX manual for details
  donor: string
    element symbol for donor ion. Currently "H" or "He" accepted.
  datacache: dict
    pass to the code to speed calculations (caches results)
  ndist : bool
    If True, calculate spectra separately for each n shell. Otherwise
    calculate entire spectrum for a given n'.
  nmax : int
    If ndist is True, this is the maximum n to include
  nmin : int
    If ndist is True, this is the minimum n to include
  debug : bool
    If true, create debugging file


  RETURNS
  -------
  none: writes a file to with the appropriate data pickled.
  """
  datacache = {}
  # get the cross section
  sigma = get_ion_sigma_acx(Z,z1, acxtype, donor, datacache=datacache,\
                            ndist = False, nmax=False, nmin=False, debug = False)

  # convert cross section to my asshole.
  sigma = no*100


def read_kronos_filedata(Z, z1, donorsym, kronosroot='.'):
  """
  Read the kronos database to get the data file name and type

  INPUTS
  ------
  Z : int
    atomic numbers
  z1 : int
    charge +1 of recombining ion
  donorsym : str
    donor atomic symbol (e.g. "H", "He")
  kronosroot : str, optional
    location of the kronos database

  RETURNS
  -------
  fname : str
    The file name (NOFILEFOUND if there is no data)
  ftype : str
    The file type: none, nl , n, depending on the resolution
  """
  foundfile=False

#  vecfac = numpy.vectorize(numpy.math.factorial)

  elsymb = pyatomdb.atomic.Ztoelsymb(Z)
  for ftype in ["rcmd","qmocc","mocc","aocc","ctmc","mclz","mclz_nres"]:
    fname = "%s/kronos/Kronos_v3/CXDatabase/Projectile_Ions/%s/Charge/%i/Targets/%s/%s%i+%s_sec_%s.cs"%\
          (kronosroot,\
           pyatomdb.atomic.Ztoelsymb(Z), z1-1, \
           donorsym,\
           pyatomdb.atomic.Ztoelsymb(Z).lower(), z1-1,\
           donorsym.lower(), ftype)

    if os.path.exists(fname):
      foundfile=True
      break
  print(fname)
  if not foundfile:
    fname = "NOFILEFOUND"
    ftype = 'None'
  else:
    if 'nres' in fname:
      ftype='n'
    else:
      ftype = 'nl'

  return fname, ftype

def calc_cross_sections(Z, z1, fname, ftype):
  """
  Calculate the cross sections into each level

  INPUTS
  ------
  Z : int
    atomic numbers
  z1 : int
    charge +1 of recombining ion
  fname : str
    filename
  kronosroot : str, optional
    filetype

  RETURNS
  -------
  fname : str
    The file name (NOFILEFOUND if there is no data)
  ftype : str
    The file type: none, nl , n, depending on the resolution
  """


  if ((ftype=='n') or (ftype=='nl')):
    dat = open(fname,'r')
    iskip = 0
    for l in dat.readlines():
      iskip+=1
      if '# (eV/u)' in l:
        lbls = l.split()[1:]
        break
    dat.close()
    d = numpy.loadtxt(fname, unpack=True, skiprows=iskip)
    energies = d[0]
    cxtype = numpy.dtype({'names':['n','l','s','CXType', 'E_min','E_max','nC', 'E', 'C'],\
                        'formats':[int, int, float,int, float, float, int, (float,d.shape[1]),\
                               (float,d.shape[1])]})
##
##
##
    if ftype=='n':
  # find all the n shells
      for i,ilbl in enumerate(lbls[1:-1]):
        #print i,ilbl ,'of', len(lbls[1:-1])

        n = int(ilbl.split('=')[1])
        l = numpy.nan
        cxtmp = numpy.zeros(1, dtype=cxtype)


        q=float(Z+1-z1)

        cxtmp['n']=n
        cxtmp['l']=-1
        cxtmp['s']=numpy.nan
        cxtmp['CXType'] = 101
        cxtmp['E_min'] = min(energies)/1000.0
        cxtmp['E_max'] = max(energies)/1000.0
        cxtmp['nC'] = len(energies)
        cxtmp['E'] = energies/1000.0
        cxtmp['C'] = d[i+1,:]*1e-16

        if i == 0:
          cxdat = cxtmp
        else:
          cxdat = numpy.append(cxdat, cxtmp)


    else:

      cxdat = numpy.zeros(d.shape[0]-2,dtype=cxtype)
      llist = 'spdfghiklmnoqrtuvwxyz'
      for i in range(len(cxdat)):
        if lbls[i+1]=='Total':
          cxdat[i]['n']=-1
          continue

        cfg,term = lbls[i+1].split('(')
        term=term[:-1]

        s=(int(term[0])-1.0)/2.0

        n= re.match('^[0-9]*',cfg).group()
        lsymb = cfg[len(n)]
        n=int(n)
        l = llist.index(lsymb)
        cxdat[i]['n'] = n
        cxdat[i]['l'] = l
        cxdat[i]['s'] = s
        cxdat[i]['CXType'] = 100
        cxdat[i]['E_min'] = min(energies)/1000.0
        cxdat[i]['E_max'] = max(energies)/1000.0
        cxdat[i]['nC'] = len(energies)
        cxdat[i]['E'] = energies/1000.0
        cxdat[i]['C'] = d[i+1,:]*1e-16
      cxdat=cxdat[cxdat['n']>0]

  return cxdat

def write_sigma_cx(sigma_cx_fname, sigma_cx, Z, z1, donorsym, ftype, fname):
  """
  Write out the charge exchange cross section file

  INPUTS
  ------
  sigma_cx_fname : string
    File name to write out to
  sigma_cx : numpy.array
    The charge exchange data from calc_cross_sections
  Z : int
    atomic numbers
  z1 : int
    charge +1 of recombining ion
  donorsymb : string
    Donor symbol
  ftype : string
    data type (nl, n, or none)
  fname : string
    file data was read from

  RETURNS
  -------
  none
  """

  hdu0 = pyatomdb.pyfits.PrimaryHDU()


  hdu1 = pyatomdb.pyfits.BinTableHDU.from_columns(pyatomdb.pyfits.ColDefs(
            [pyatomdb.pyfits.Column(name='n',
               format='1J',
               array=sigma_cx['n']),
             pyatomdb.pyfits.Column(name='l',
               format='1J',
               array=sigma_cx['l']),
             pyatomdb.pyfits.Column(name='s',
               format='1E',
               array=sigma_cx['s']),
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
  if donorsym == 'H':
    donormass = 1.0079
  elif donorsym == 'He':
    donormass = 4.0026
  elif donorsym == 'CO':
    donormass = 28.0101
  elif donorsym == 'CO2':
    donormass = 44.0095
  elif donorsym == 'H2':
    donormass = 2.01589
  elif donorsym == 'H20':
    donormass = 18.0153
  elif donorsym == 'N2':
    donormass = 28.0134
  else:
    print("UNKNOWN DONOR SYSTEM: %s"%(donor))
    exit

  hdu1.header['DonMass']=(donormass, "Donor Mass in a.m.u.")
  hdu1.header['RecMass']=(pyatomdb.atomic.Z_to_mass(Z), "Receiver Mass in a.m.u.")
  hdu1.header['System']="%s+%i + %s"%(pyatomdb.atomic.Ztoelsymb(Z), z1-1, donorsym)

  hdu1.header['Source']='Kronos_v3'
  tmpfname = fname
  while len(tmpfname)>40:
    hdu1.header.add_comment('SRC FILE: %s'%(tmpfname[:40]))
    tmpfname = tmpfname[40:]
  hdu1.header.add_comment('SRC FILE: %s'%(tmpfname))

  hdulist = pyatomdb.pyfits.HDUList([hdu0,hdu1])

  s = '.'
  ss = sigma_cx_fname.split('/')
  if len(ss) > 1:

    for sss in ss[:-1]:
      s+='/%s'%(sss)
      pyatomdb.util.mkdir_p(s)

  hdulist.writeto(sigma_cx_fname, checksum=True, overwrite=True)

def generate_datatypes(dtype, npseudo=0, ncontinuum=0):
  """
  returns the various data types needed by apec

  Parameters
  ----------
  dtype : string
    One of "linetype", "cielinetype", "continuum"
  npseudo : int (default=0)
    Number of pseudocontinuum points for "continuum" type
  ncontinuum : int (default=0)
    Number of continuum points for "continuum" type


  Returns
  -------
  numpy.dtype
    The data dtype in question
  """
  if dtype == 'linetype':
    ret = numpy.dtype({'names':['lambda',\
                                 'lambda_err',\
                                 'epsilon',\
                                 'epsilon_err',\
                                 'element',\
                                 'ion', \
                                 'elem_drv',\
                                 'ion_drv', \
                                 'upperlev',\
                                 'lowerlev'],\
                        'formats':[float,\
                                   float,\
                                   float,\
                                   float,\
                                   int,\
                                   int,\
                                   int,\
                                   int,\
                                   int,\
                                   int]})

  elif dtype == 'linetype_cx_nl':
    ret = numpy.dtype({'names':['lambda',\
                                'epsilon',\
                                'element',\
                                'ion', \
                                'upperlev',\
                                'lowerlev',\
                                'n',\
                                'l'],\
                       'formats':[float,\
                                  float,\
                                  int,\
                                  int,\
                                  int,\
                                  int,\
                                  int,\
                                  int]})

  elif dtype =='linelist_cie':
    ret = numpy.dtype({'names':['lambda',\
                                 'lambda_err',\
                                 'epsilon',\
                                 'epsilon_err',\
                                 'element',\
                                 'ion', \
                                 'upperlev',\
                                 'lowerlev'],\
                        'formats':[float,\
                                   float,\
                                   float,\
                                   float,\
                                   int,\
                                   int,\
                                   int,\
                                   int]})

  elif dtype == 'linetype_cap':
    ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'Elem_Drv',\
                                 'Ion_Drv', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[float,\
                                   float,\
                                   float,\
                                   float,\
                                   int,\
                                   int,\
                                   int,\
                                   int,\
                                   int,\
                                   int]})

  elif dtype == 'linelist_cie_cap':
    ret = numpy.dtype({'names':['Lambda',\
                                 'Lambda_Err',\
                                 'Epsilon',\
                                 'Epsilon_Err',\
                                 'Element',\
                                 'Ion', \
                                 'UpperLev',\
                                 'LowerLev'],\
                        'formats':[float,\
                                   float,\
                                   float,\
                                   float,\
                                   int,\
                                   int,\
                                   int,\
                                   int]})


  elif dtype == 'continuum':
    if ncontinuum==0:
      ncontinuum+=1
    if npseudo==0:
      npseudo+=1

    ret = numpy.dtype({'names':['Z','rmJ','N_Cont','E_Cont','Continuum','Cont_Err','N_Pseudo','E_Pseudo','Pseudo','Pseudo_Err'],\
                       'formats':[int, int, \
                                  int, (float, ncontinuum), (float, ncontinuum),(float, ncontinuum),\
                                  int, (float, npseudo), (float, npseudo),(float, npseudo)]})
  elif dtype == 'lineparams':
    ret = numpy.dtype({'names':['kT','EDensity','Time','Nelement','Nline'],\
                       'formats':[float, float, float, int, int]})
  elif dtype == 'cocoparams':
    ret = numpy.dtype({'names':['kT','EDensity','Time','NElement','NCont', 'NPseudo'],\
                       'formats':[float, float, float, int, int, int]})
  elif dtype == 'ecdat':
    # Electron collisional data
    ret = numpy.dtype({'names':['lower_lev','upper_lev','coeff_type',\
                                'min_temp','max_temp', 'temperature', \
                                'effcollstrpar','inf_limit','reference'],\
                       'formats':[int, int, int, \
                                  float, float, (float,20), \
                                  (float,20), float, '|S20']})


  else:
    print("Unknown dtype %s in generate_datatypes"%(dtype))
  return ret


def gen_configs(lvdat):
  # this routine gets a list of all the unique configurations in the file
  cfglist = pyatomdb.util.unique(lvdat[1].data.field('elec_config'))

  return cfglist

def check_1_config_match(cfg1, cfg2):
  # this routine matches 2 occupancy vectors, true if match, false otherwise
  #print("cfg1:", cfg1, ", cfg2:", cfg2)
  i1 = numpy.where(cfg1 > 0)[0]
  match=True
  for i in i1:
    if cfg1[i] != cfg2[i]:
      match=False
  if sum(cfg1) != sum(cfg2):
    match=False

  return match

def check_config_match(cfg, cfglist):
  # this routine looks for a matching occupancy vector in a supplied set
  match=False
  imatch = -1
  for icfg in range(len(cfglist)):
    if check_1_config_match(cfg, cfglist[icfg]):
      match = True
      imatch=icfg
      break
  return match, imatch

def gen_autos_input(lvconfigs, gndconfig, nmax, nel, Z, z1r):
  """
   this routine will generate an AUTOSTRUCTURE input file containing
   all the configurations in the lvconfigs array AND all the singly
   excited configs up to and including nmax

  PARAMETERS
  ----------
  lvconfigs: list[strings]
    The list of configurations in the APEC file
    (e.g. ['1s2 2p1', '1s2 2s1'])
  gndconfig: string
    The ground configuration (e.g. '1s2 2s1')
  nmax : int
    The maximum n shell to promote to.
  nel : int
    Number of electrons in RECOMBINED ion
  Z : int
    Nuclear charge
  z1r : int
    Ion charge +1 of the RECOMBINED ion

  RETURNS
  -------
  dirname : string
    The directory where the input file has been written
  """

  import shutil

  cfg = []

  shlmax=sum(range(1,nmax+1))
  maxoccup = numpy.zeros(shlmax)
  maxoccup[0] = 2
  ntmp = 1
  ltmp = 0
  for i in range(1,shlmax):
    if ltmp == ntmp-1:
      ntmp +=1
      ltmp = 0
    else:
      ltmp += 1
    maxoccup[i] = ltmp*4+2

  groundconfig, zzz = pyatomdb.atomic.config_to_occup(gndconfig, nel, shlmax=shlmax)
  occlist=[]
  occlist.append(groundconfig)
  for i in range(len(lvconfigs)):
    tmpocc,zzz = pyatomdb.atomic.config_to_occup(lvconfigs[i], nel, shlmax=shlmax)

    match, imatch  = check_config_match(tmpocc, occlist)
    if match == False:
      occlist.append(tmpocc)
  # add the singly excited levels
  ival = max(numpy.where(groundconfig > 0)[0])

  if nel==1:
    groundconfig = numpy.zeros(shlmax, dtype=int)
  else:
    plvdat  = pyatomdb.atomdb.get_data(Z, z1r+1, 'LV')
    gndcfg = plvdat[1].data['ELEC_CONFIG'][0]
    groundconfig,zzz = pyatomdb.atomic.config_to_occup(gndcfg, nel-1, shlmax=shlmax)



  for i in range(shlmax):
    tmpocc = groundconfig.copy()
    if tmpocc[i] == maxoccup[i]: continue

    tmpocc[i]+=1
    match, imatch  = check_config_match(tmpocc, occlist)
    if match == False:
      occlist.append(tmpocc)

  mxconf = len(occlist)
  mxvorb = shlmax
  elsymb = pyatomdb.atomic.z0toelsymb(Z)
  dirname ='autos/'+elsymb.lower()+repr(z1r)

  # remove anything hanging around in the autostructure directory
  shutil.rmtree(dirname, ignore_errors=True)

  pyatomdb.util.mkdir_p('autos')
  pyatomdb.util.mkdir_p(dirname)
  print("HIINPUTWRITES")
  fname = dirname+'/input.dat'
  f=open(fname, 'w')
  f.write('A.S. '+elsymb+'+'+repr(z1-1)+' structure\n')
  f.write(' &SALGEB MXCONF='+repr(mxconf)+' MXVORB='+repr(mxvorb)+\
          ' RAD=\'E2\' CUP=\'ICR\' &END\n')

  s_out = ''
  l = 0
  n = 1
  for i in range(shlmax):
    s_out=s_out+"%2i %2i " %(n,l)
    if n-l==1:
      n = n+1
      l = 0
    else:
      l = l +1
  s_out=s_out[:-1]+'\n'
  f.write(s_out)

  for iocc in occlist:
    f.write(occ_to_autos_string(iocc))

  f.write(' &SMINIM NZION='+repr(Z)+ ' &END\n')

  f.close()


  return dirname

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def occ_to_autos_string(iocc):
  """
  Convert occupancy (in array of integers, from 1s, 2s, 2p, 3s etc) to
  an Autostructure input array.

  PARAMETERS
  ----------
  iocc : array[int]
    Occupancies, e.g. [2,0,1,0,1] for 1s2 2p1 3p1

  RETURNS
  -------
  s : string
    The occupancy as autostructure string "  2   0   1   0   1"
  """

  s=''
  for i in iocc:
    s = s + " %3i  " % (i)
  s=s[:-2]+'\n'
  return s

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def generate_autostructure(Z, z1r, nmax, lvdat=False, filemap = False, outdir='autos'):
  """
  Generate the autostructure input files and run the code

  PARAMETERS
  ----------
  Z : int
    nuclear charge
  z1r : int
    ion charge +1 of the *recombined* ion
  lvdat : hdulist
    the energy level file already opened
  nmax : int
    max n shell to go to
  filemap : filemap
    the filemap to update
  outdir : directory (optional)
    Autostructure files will be in outdir/ELSYMB/ELSYMB_Z1

  RETURNS
  -------
  None
  """
  import autos, shutil

# get the LV configs
#  print lvfile
  if lvdat != False:

    nmax2 = numpy.max(lvdat[1].data.field('n_quan'))

    gndconfig = lvdat[1].data.field('elec_config')[0]
    lvconfigs = gen_configs(lvdat)


    ### HERE ###
  else:
    nmax2 = 0
     # get ground config from adas
    a00file = '../adf00/'+pyatomdb.atomic.z0toelsymb(Z).lower()+'.dat'
    gndconfig=pyatomdb.atomic.fmtconf(adas.read_adf00(a00file)[z1r]['gndcfg'])
    lvconfigs=[gndconfig]

  nel = Z+1-z1r
  ndir = gen_autos_input(lvconfigs, gndconfig, max([nmax,nmax2]), nel, Z, z1r)

  # change directory, run autostructure?
  ##HYDRACHANGE
  #hydra
  autosbin = '/home/afoster/autostructure/autos.x'
  #laptop
  #autosbin='/export1/autostructure/autos.x'
  autosbin='/home/afoster/work/export1/projects/dr/codes/aslm.x'
  odir = os.getcwd()
  #ndir = '../autos/'+pyatomdb.atomic.z0toelsymb(z0).lower()+repr(z1)


  # make autostructure directory
  pyatomdb.util.mkdir_p(ndir)

  # change into it
  print("changing into %s"%(ndir))
  os.chdir(ndir)
  print("now in %s"%( os.getcwd()))

  # run autostructure
  cmd = autosbin+' < input.dat'
  print("Running: %s in %s"%(cmd, ndir))
  os.system(cmd)

  # change back to initial directory
  os.chdir(odir)


  # convert autostructure to FITS files

  autdat = autos.read_oic(ndir)
  #print(autdat.keys(0))
  autdat = autdat[list(autdat.keys())[0]]
  autdat = autdat[list(autdat.keys())[0]]


  # create temporary FITS file names
  tmplvfilename = ndir+'/%s_%i_LV_autos.fits'%(pyatomdb.atomic.Ztoelsymb(Z).lower(), z1r)
  tmplafilename = ndir+'/%s_%i_LA_autos.fits'%(pyatomdb.atomic.Ztoelsymb(Z).lower(), z1r)

  # assemble data

  lvtmp = {}
  lvtmp['Z'] = Z
  lvtmp['z1'] = z1r
  lvtmp['comments'] = ['Generated by Autostructure for ACX 2 run']
  #for ik, k in enumerate(autdat.keys()):
    #print("%i : "%(ik), k)
  lv_to_k = numpy.zeros(len(autdat['lev'])+1, dtype=int)
  for i in range(len(autdat['lev'])):
    lv_to_k[autdat['lev']['lv'][i]] = autdat['lev']['k'][i]


  lvdat = numpy.zeros(len(autdat['lev']), dtype = numpy.dtype({\
                      'names':['elec_config',\
                               'energy',\
                              'e_error',\
                              'n_quan',\
                              'l_quan',\
                              's_quan',\
                              'lev_deg',\
                              'phot_type',\
                              'phot_par',\
                              'aaut_tot',\
                              'arad_tot',\
                              'energy_ref',\
                              'phot_ref',\
                              'aaut_ref',\
                              'arad_ref'],\
                      'formats': ['|S40',\
                                  float,\
                                  float,\
                                  int,\
                                  int,\
                                  float,\
                                  int,\
                                  int,\
                                  (float, 20),\
                                  float,\
                                  float,\
                                  '|S20',\
                                  '|S20',\
                                  '|S20',\
                                  '|S20']}))

  lvdat['elec_config'] = autdat['lev']['cfgstr']
  lvdat['energy'] = autdat['lev']['energy']
  lvdat['e_error'][:] = numpy.nan


  for i in range(len(lvdat)):
    lvdat['n_quan'][i] = pyatomdb.atomdb._extract_n(lvdat['elec_config'][i].decode('ascii'))

  lvdat['l_quan'] = autdat['lev']['l']
  lvdat['s_quan'] = (autdat['lev']['s2p1']-1.0)/2.0
  lvdat['lev_deg'] = autdat['lev']['j2']+1
  lvdat['phot_type'] = -1

  for i in range(len(lvdat)):
    ii = numpy.where(autdat['trn']['upper_lev']==autdat['lev']['lv'][i])[0]

    if len(ii) > 0:
      lvdat['arad_tot'][i] = sum(autdat['trn']['ar'][ii])
  lvdat['energy_ref'][:] = 'ACX_Autostr'
  lvdat['arad_ref'][:] = 'ACX_Autostr'

  lvtmp['data'] = lvdat

  pyatomdb.util.write_lv_file(tmplvfilename, lvtmp, clobber=True)

  # THIS IS BAWSAX.
  latmp = {}
  latmp['Z'] = Z
  latmp['z1'] = z1r
  latmp['comments'] = ['Generated by Autostructure for ACX 2 run']

  ladat=numpy.zeros(len(autdat['trn']), dtype = numpy.dtype({\
                        'names': ['upper_lev',\
                                  'lower_lev',\
                                  'wavelen',\
                                  'wave_obs',\
                                  'wave_err',\
                                  'einstein_a',\
                                  'ein_a_err',\
                                  'wave_ref',\
                                  'wv_obs_ref',\
                                  'ein_a_ref'],\
                        'formats' :[int,\
                                    int,\
                                    float,\
                                    float,\
                                    float,\
                                    float,\
                                    float,\
                                    '|S20',\
                                    '|S20',\
                                    '|S20']}))
  ladat['upper_lev'] = lv_to_k[autdat['trn']['upper_lev']]
  ladat['lower_lev'] = lv_to_k[autdat['trn']['lower_lev']]

  ladat['wavelen'] = pyatomdb.const.HC_IN_KEV_A/\
                     (autdat['trn']['de_ryd']*pyatomdb.const.RYDBERG)

  ladat['wave_err'][:] = numpy.nan
  ladat['wave_obs'][:] = numpy.nan
  ladat['einstein_a'] = autdat['trn']['ar']
  ladat['ein_a_err'][:] = numpy.nan
  ladat['wave_ref'][:] = 'ACX_Autostr'
  ladat['ein_a_ref'][:] = 'ACX_Autostr'

  latmp['data'] = ladat

  pyatomdb.util.write_la_file(tmplafilename, latmp, clobber=True)

  return tmplvfilename, tmplafilename


def get_n(cfg):
  icfg = cfg.split()
  n=-1
  for i in icfg:
    ntmp = int(re.match('[0-9]+',i).group(0))
    n = max([n,ntmp])
  return n

def fix_jj_doubles(cfgstr):
  max_n = get_n(cfgstr)
  occup = numpy.zeros( sum(range(1,max_n+1)), dtype=int)
#  print occup
  llist='spdfghiklmnoqrtuvwxyz'

  cfg = cfgstr.split(' ')
  for icfg in cfg:
    ntmp = re.search("^[0-9]+",icfg)
    n=int(ntmp.group(0))
    ltmp = re.search("[a-zA-Z]",icfg)
    l=llist.index(ltmp.group(0))
    otmp = re.search("[0-9]+$",icfg)
    o=int(otmp.group(0))

    # find shell for this nl combo
    shlind = sum(range(n))+l
    occup[shlind] += o

  # convert back to string
  s = ''
  n = 1
  l = 0
  for i in range(len(occup)):
    if occup[i] != 0:
      s = s +repr(n)+llist[l]+repr(occup[i])+' '
    if n-l==1:
      n = n + 1
      l = 0
    else:
      l = l + 1
  s = s[:-1]
  return s

def check_atomic_data(Z, z1, max_n, fmapfile=False):
  """
  Check that we have sufficiently high n-shell atomic data

  PARAMETERS
  ----------
  Z: int
    atomic number
  z1 : int
    nuclear charge of recombining ion +1
  max_n : int
    the maximum n shell we require data for
  fmapfile : string (optional)
    The name of the file. If set, will also be updated with newly
    generated file
  """

  # read the atomic data from existing APEC files.
  fmap = pyatomdb.atomdb.read_filemap(fmapfile)

  # find the corresponding energy level file
  # since z1 is recombining, we are looking for z1-1
  ifmap = numpy.where((fmap['Z']==Z) & (fmap['z1']==z1-1))[0][0]
  lvfile = fmap['lv'][ifmap]
  lafile = fmap['la'][ifmap]

  try:
    lvfile = lvfile.decode('ascii')
    lafile = lafile.decode('ascii')
  except:
    pass

  finallvfilename = ''
  finallafilename = ''
  # check if there is a file. If so, check the maximum N shell.
  if not lvfile=='':
    #help(lvfile)
    lvdat = pyatomdb.pyfits.open(lvfile)

#    if (((Z==23)|(Z==21)) & (z1==7)):
#      print("Hack for v7")
#      zzz=input('argh')
#      autlvfilename, autlafilename= generate_autostructure(Z, z1-1, max_n, lvdat=lvdat, filemap=fmapfile)

#      finallvfilename, finallafilename = \
#         merge_autostructure(Z, z1-1, autlvfilename, autlafilename, lvfile, lafile, filemap=fmapfile )


    if max(lvdat[1].data['n_quan']) >= max_n:
      # all is well
      pass

    else:
      print("Insufficient n-shell data available in %s (need %i, have %i). Generating atomic data for Z=%i, z1=%i"%\
           (lvfile, max_n, max(lvdat[1].data['n_quan']),Z, z1-1))

      autlvfilename, autlafilename= generate_autostructure(Z, z1-1, max_n, lvdat=lvdat, filemap=fmapfile)

      finallvfilename, finallafilename = \
         merge_autostructure(Z, z1-1, autlvfilename, autlafilename, lvfile, lafile, filemap=fmapfile )

  else:
    print("No APED data available for Z=%i, z1=%i, running Autostructure."%(Z, z1-1))

    finallvfilename, finallafilename= generate_autostructure(Z, z1-1, max_n, lvdat=False, filemap=fmapfile)

  if finallvfilename != '':

    fmap = pyatomdb.atomdb.read_filemap(fmapfile)

  # find the corresponding energy level file
  # since z1 is recombining, we are looking for z1-1
    print("finallvfilename=%s"%( finallvfilename))
    ifmap = numpy.where((fmap['Z']==Z) & (fmap['z1']==z1-1))[0][0]
    fmap['lv'][ifmap] = finallvfilename
    fmap['la'][ifmap] = finallafilename

    pyatomdb.atomdb.write_filemap(fmap, fmapfile)




def merge_autostructure(Z, z1r, autlvfilename, autlafilename, lvfname, lafname, filemap=False):
  """
  Merge the autostructure results into the APED database

  PARAMETERS
  ----------
  Z : int
    atomic number
  z1r : int
    nuclear charge of recombined ion +1
  autlvfilename : string
    file name of Autostructure LV data
  autlafilename : string
    file name of Autostructure LA data
  fmapfile : string (optional)
    The name of the file. If set, will also be updated with newly
    generated file
  """

  apedlvdat = pyatomdb.pyfits.open(lvfname)
  autlvdat = pyatomdb.pyfits.open(autlvfilename)

  apedladat = pyatomdb.pyfits.open(lafname)
  autladat = pyatomdb.pyfits.open(autlafilename)

  # match the levels
  nel=Z+1-z1r
  if sum(apedlvdat[1].data['l_quan']<0) > 0:
    jjcup=True
  else:
    jjcup = False

  for i in range(len(apedlvdat[1].data)):
    tmp,x = pyatomdb.atomic.config_to_occup(apedlvdat[1].data['elec_config'][i],\
                                            nel=nel)


    apedlvdat[1].data['elec_config'][i] = pyatomdb.atomic.occup_to_cfg(tmp)

  if jjcup:
    for i in range(len(apedlvdat[1].data)):

      tmp = fix_jj_doubles(apedlvdat[1].data['elec_config'][i])
      tmp,x = pyatomdb.atomic.config_to_occup(tmp,\
                                            nel=nel)

      apedlvdat[1].data['elec_config'][i] = pyatomdb.atomic.occup_to_cfg(tmp)

  auttolv = numpy.zeros(len(autlvdat[1].data), dtype=int)
  lvtoaut = numpy.zeros(len(apedlvdat[1].data), dtype=int)

  auttolv[:]=-1
  lvtoaut[:]=-1
  badcfg = []

  for ilv in range(len(lvtoaut)):
    if lvtoaut[ilv] != -1: continue
    if ((apedlvdat[1].data['l_quan'][ilv]<0) | (apedlvdat[1].data['s_quan'][ilv]<0)):
      levjjcup=True
    else:
      levjjcup=False

    if levjjcup:
      i = numpy.where((apedlvdat[1].data['elec_config']==\
                           apedlvdat[1].data['elec_config'][ilv]) &\
                      (apedlvdat[1].data['lev_deg']==\
                           apedlvdat[1].data['lev_deg'][ilv]))[0]

      j = numpy.where((autlvdat[1].data['elec_config']== \
                           apedlvdat[1].data['elec_config'][ilv]) &\
                      (autlvdat[1].data['lev_deg']==\
                           apedlvdat[1].data['lev_deg'][ilv]))[0]
      print('is jj')

    else:
      i = numpy.where((apedlvdat[1].data['elec_config']==\
                           apedlvdat[1].data['elec_config'][ilv]) &\
                      (apedlvdat[1].data['l_quan']==\
                           apedlvdat[1].data['l_quan'][ilv]) &\
                      (apedlvdat[1].data['s_quan']==\
                           apedlvdat[1].data['s_quan'][ilv]) &\
                      (apedlvdat[1].data['lev_deg']==\
                           apedlvdat[1].data['lev_deg'][ilv]))[0]

      j = numpy.where((autlvdat[1].data['elec_config']== \
                           apedlvdat[1].data['elec_config'][ilv]) &\
                      (autlvdat[1].data['l_quan']== \
                           apedlvdat[1].data['l_quan'][ilv]) &\
                      (autlvdat[1].data['s_quan']==\
                           apedlvdat[1].data['s_quan'][ilv]) &\
                      (autlvdat[1].data['lev_deg']==\
                           apedlvdat[1].data['lev_deg'][ilv]))[0]
      print('not jj')

    print('found %i matches in APED, %i in AUTOS'%(len(i), len(j)))
    for iii in i:
      print("%6i %40s %e %3i %.1f %3i"%\
                  (iii+1,
                   apedlvdat[1].data['elec_config'][iii], \
                   apedlvdat[1].data['energy'][iii],\
                   apedlvdat[1].data['l_quan'][iii],\
                   apedlvdat[1].data['s_quan'][iii],\
                   apedlvdat[1].data['lev_deg'][iii]))
    print('--')
    for iii in j:
      print("%6i %40s %e %3i %.1f %3i"%\
                  (iii+1,
                   autlvdat[1].data['elec_config'][iii], \
                   autlvdat[1].data['energy'][iii],\
                   autlvdat[1].data['l_quan'][iii],\
                   autlvdat[1].data['s_quan'][iii],\
                   autlvdat[1].data['lev_deg'][iii]))
    print('XX')

    if (len(j) != len(i)):
      alreadydone = False
      for ibad in badcfg:
        if ilv in ibad['aped']:
          alreadydone=True
      if alreadydone: continue
      print('WARNING: mismatch in number of matches for :')
      print("%20s %2i %3.1f %2i" %(apedlvdat[1].data['elec_config'][ilv],\
                                   apedlvdat[1].data['l_quan'][ilv],\
                                   apedlvdat[1].data['s_quan'][ilv],\
                                   apedlvdat[1].data['lev_deg'][ilv]))

      print(" APED: %i, AUTOS: %i" %(len(i), len(j)))
      print("  APEDind: ",i)
      print("  AUTOSind: ",j)
      badcfg.append({})
      badcfg[-1]['elec_config'] = apedlvdat[1].data['elec_config'][ilv]
      badcfg[-1]['l_quan'] = apedlvdat[1].data['l_quan'][ilv]
      badcfg[-1]['s_quan'] = apedlvdat[1].data['s_quan'][ilv]
      badcfg[-1]['lev_deg'] = apedlvdat[1].data['lev_deg'][ilv]
      badcfg[-1]['aped']=i
      badcfg[-1]['autos']=j
      badcfg[-1]['parity'] = \
               pyatomdb.atomic.get_parity(apedlvdat[1].data['elec_config'][ilv])
      badcfg[-1]['done'] = False


    elif (len(j) > 1):
      ii = i[numpy.argsort(apedlvdat[1].data['energy'][i])]
      jj = j[numpy.argsort(autlvdat[1].data['energy'][j])]

      for k in range(len(ii)):
        auttolv[jj]=ii
        lvtoaut[ii]=jj

    elif (len(j) == 1):
      auttolv[j[0]]=i[0]
      lvtoaut[i[0]]=j[0]

    else:
      auttolv[j[0]]=i[0]
      lvtoaut[i]=-2
#    print auttolv
#    print lvtoaut
    # deal with "badcfg" - those with different numbers of levels
  if len(badcfg) > 0:
#      print badcfg
    for ibad, bad in enumerate(badcfg):
      if ibad==len(badcfg)-1:
        continue
      if bad['done']== True: continue
      badcfg[ibad]['done']= True
      matches= [ibad]
      for iibad in range(ibad+1, len(badcfg)):
        if badcfg[iibad]['done']== True: continue
#          if (jjcup):
        if ((badcfg[iibad]['lev_deg'] ==bad['lev_deg']) &\
            (badcfg[iibad]['parity'] ==bad['parity'])):
          matches.append(iibad)
          badcfg[iibad]['done']= True
#          else:
#            if ((badcfg[iibad]['lev_deg'] ==bad['lev_deg']) &\
#                (badcfg[iibad]['l_quan'] ==bad['l_quan']) &\
#                (badcfg[iibad]['s_quan'] ==bad['s_quan']) &\
#                (badcfg[iibad]['parity'] ==bad['parity'])):
#              matches.append(iibad)
#              print "adding match: %i"%(iibad)
#              badcfg[iibad]['done']= True
      matches_aped=[]
      matches_autos=[]
#        print "matches:",matches
      for m in matches:
#          print "m=",m
        badcfg[m]['done']=True
#          print "badcfg[%i]['aped']="%(m), badcfg[m]['aped']
        for mm in badcfg[m]['aped']:

          matches_aped.append(mm)
        for mm in badcfg[m]['autos']:
          matches_autos.append(mm)
      matches_aped.sort()
      matches_autos.sort()

      print("matches_aped:", matches_aped  )
      print("matches_autos:", matches_autos)

#        matches_aped = bad['aped']
#        matches_autos = bad['autos']

      for ii in range(max([len(matches_aped), len(matches_autos)])):
        if ((ii < len(matches_aped)) & (ii < len(matches_autos))):
          lvtoaut[matches_aped[ii]] = matches_autos[ii]
          auttolv[matches_autos[ii]] = matches_aped[ii]
        elif ((ii < len(matches_aped)) & (ii >= len(matches_autos))):
          print("ERROR: no match in AUTOS found for APED level %i"%\
                 (matches_aped[ii]))
        elif ((ii >= len(matches_aped)) & (ii < len(matches_autos))):
          #lvtoaut[matches_aped[ii]] = matches_autos[ii]
          auttolv[matches_autos[ii]] = -2

  n_levels = len(lvtoaut)+len(numpy.where(auttolv<0)[0])

  for i in range(len(lvtoaut)):
    print("lvtoaut[%i] = %i"%(i+1, lvtoaut[i]+1))

  # dump evel translations to file
  d = {}
  d['lvtoaut'] = lvtoaut
  d['auttolv'] = auttolv

#  pickle.dump(d, open('pickles/dumplvtrans_%i_%i.pkl'%(Z,z1), 'wb'))
#  exit()

  lvdat=numpy.zeros(n_levels, dtype=
                    numpy.dtype({'names':['elec_config','energy',\
                                          'e_error', 'n_quan',\
                                          'l_quan', 's_quan',\
                                          'lev_deg', 'phot_type', \
                                          'phot_par', 'aaut_tot',\
                                          'arad_tot','energy_ref', \
                                          'phot_ref', 'aaut_ref',\
                                          'arad_ref'],\
                                 'formats':['S40', float,\
                                            float, int, \
                                            int, float,\
                                            int, int, \
                                            (float,20), float,\
                                            float, '|S20',\
                                            '|S20','|S20',\
                                            '|S20']}))

  lvdat2={}
  lvdat2['Z']=Z
  lvdat2['z1']=z1r
  lvdat2['comments']=['made for CX studies only']

  for i in range(len(lvtoaut)):
    lvdat['elec_config'][i] = apedlvdat[1].data['elec_config'][i]
    lvdat['energy'][i]      = apedlvdat[1].data['energy'][i]
    lvdat['e_error'][i]     = apedlvdat[1].data['e_error'][i]
    lvdat['n_quan'][i]      = apedlvdat[1].data['n_quan'][i]
    lvdat['l_quan'][i]      = apedlvdat[1].data['l_quan'][i]
    lvdat['s_quan'][i]      = apedlvdat[1].data['s_quan'][i]
    lvdat['lev_deg'][i]     = apedlvdat[1].data['lev_deg'][i]
    lvdat['phot_type'][i]   = apedlvdat[1].data['phot_type'][i]
    lvdat['phot_par'][i]    = apedlvdat[1].data['phot_par'][i]
    lvdat['energy_ref'][i]  = apedlvdat[1].data['energy_ref'][i]
    lvdat['phot_ref'][i]    = apedlvdat[1].data['phot_ref'][i]
    if 'ARAD_TOT' in apedlvdat[1].data.names:
      lvdat['aaut_tot'][i]   = apedlvdat[1].data['aaut_tot'][i]
      lvdat['arad_tot'][i]    = apedlvdat[1].data['arad_tot'][i]
      lvdat['arad_ref'][i]  = apedlvdat[1].data['arad_ref'][i]
      lvdat['aaut_ref'][i]    = apedlvdat[1].data['aaut_ref'][i]
    else:
      lvdat['arad_ref'][i]  = 'FIXME'
      lvdat['aaut_ref'][i]    = 'FIXME'

  inext = len(lvtoaut)

  j = numpy.where(auttolv < 0)[0]

  for i in j:
    lvdat['elec_config'][inext] = autlvdat[1].data['elec_config'][i]
    lvdat['energy'][inext]      = autlvdat[1].data['energy'][i]
    lvdat['e_error'][inext]     = numpy.nan
    lvdat['n_quan'][inext]      = autlvdat[1].data['n_quan'][i]
    lvdat['l_quan'][inext]      = autlvdat[1].data['l_quan'][i]
    lvdat['s_quan'][inext]      = autlvdat[1].data['s_quan'][i]
    lvdat['lev_deg'][inext]     = autlvdat[1].data['lev_deg'][i]
    lvdat['phot_type'][inext]   = autlvdat[1].data['phot_type'][i]
    lvdat['phot_par'][inext]    = autlvdat[1].data['phot_par'][i]
    lvdat['energy_ref'][inext]  = autlvdat[1].data['energy_ref'][i]
    lvdat['phot_ref'][inext]    = autlvdat[1].data['phot_ref'][i]
    lvdat['aaut_tot'][inext]   = autlvdat[1].data['aaut_tot'][i]
    lvdat['arad_tot'][inext]    = autlvdat[1].data['arad_tot'][i]
    lvdat['arad_ref'][inext]  = autlvdat[1].data['arad_ref'][i]
    lvdat['aaut_ref'][inext]    = autlvdat[1].data['aaut_ref'][i]
    auttolv[i] = inext
    inext = inext + 1

  lvdat2['data'] = numpy.array(lvdat)


  # Now do the LA file

  ntran_max = len(apedladat[1].data) + len(autladat[1].data)
  ladat=numpy.zeros(ntran_max, dtype=\
                    numpy.dtype({'names':['upper_lev', 'lower_lev', \
                                          'wavelen', 'wave_obs', \
                                          'wave_err', 'einstein_a', \
                                          'ein_a_err', 'wave_ref', \
                                          'wv_obs_ref', 'ein_a_ref'],\
                                 'formats':[int, int,\
                                            float, float,\
                                            float, float,\
                                            float, '|S20',\
                                            '|S20','|S20']}))

  ladat2={}
  ladat2['Z']=Z
  ladat2['z1']=z1r
  ladat2['comments']=['made for CX studies only']

  for i in range(len(apedladat[1].data)):
    ladat['upper_lev'][i] = apedladat[1].data['upper_lev'][i]
    ladat['lower_lev'][i] = apedladat[1].data['lower_lev'][i]
    ladat['wavelen'][i]   = apedladat[1].data['wavelen'][i]
    ladat['wave_obs'][i]  = apedladat[1].data['wave_obs'][i]
    ladat['wave_err'][i]  = apedladat[1].data['wave_err'][i]
    ladat['einstein_a'][i]= apedladat[1].data['einstein_a'][i]
    ladat['ein_a_err'][i] = apedladat[1].data['ein_a_err'][i]
    ladat['wave_ref'][i]  = apedladat[1].data['wave_ref'][i]
    ladat['wv_obs_ref'][i]= apedladat[1].data['wv_obs_ref'][i]
    ladat['ein_a_ref'][i] = apedladat[1].data['ein_a_ref'][i]


  inext = len(apedladat[1].data)

  rtrckfname = 'runtrack_%i_%i.log'%(Z,z1r)
  rtrck = open(rtrckfname, 'w')
  rtrck.write('starting\n')
  rtrck.close()


  # change the levels in the aut file

  autladat[1].data['upper_lev'] -= 1
  autladat[1].data['lower_lev'] -= 1

  autladat[1].data['upper_lev'] = auttolv[autladat[1].data['upper_lev']]+1
  autladat[1].data['lower_lev'] = auttolv[autladat[1].data['lower_lev']]+1

  copyme = numpy.ones(len(autladat[1].data), dtype=bool)

  uplist = pyatomdb.util.unique(apedladat[1].data['upper_lev'])
  for up in uplist:
    lolist = pyatomdb.util.unique(apedladat[1].data['lower_lev'][apedladat[1].data['upper_lev']==up])
    iiup = numpy.where(autladat[1].data['upper_lev']==up)[0]
    if len(iiup)>0:
      for ii in iiup:
        if autladat[1].data['lower_lev'][ii] in lolist:
          copyme[ii] = False


  uplist = pyatomdb.util.unique(apedladat[1].data['lower_lev'])

  # cancel reverse transactions
  for up in uplist:
    lolist = pyatomdb.util.unique(apedladat[1].data['upper_lev'][apedladat[1].data['lower_lev']==up])
    iiup = numpy.where(autladat[1].data['upper_lev']==up)[0]
    if len(iiup)>0:
      for ii in iiup:
        if autladat[1].data['lower_lev'][ii] in lolist:
          copyme[ii] = False


  print("Including %i of %i transitions from CX calculation"%(sum(copyme), len(copyme)))
  ifin = inext+sum(copyme)

  ladat['upper_lev'][inext:ifin] = autladat[1].data['upper_lev'][copyme]
  ladat['lower_lev'][inext:ifin] = autladat[1].data['lower_lev'][copyme]
  ladat['wavelen'][inext:ifin]   = autladat[1].data['wavelen'][copyme]
  ladat['wave_obs'][inext:ifin]  = autladat[1].data['wave_obs'][copyme]
  ladat['wave_err'][inext:ifin]  = autladat[1].data['wave_err'][copyme]
  ladat['einstein_a'][inext:ifin]= autladat[1].data['einstein_a'][copyme]
  ladat['ein_a_err'][inext:ifin] = autladat[1].data['ein_a_err'][copyme]
  ladat['wave_ref'][inext:ifin]  = autladat[1].data['wave_ref'][copyme]
  ladat['wv_obs_ref'][inext:ifin]= autladat[1].data['wv_obs_ref'][copyme]
  ladat['ein_a_ref'][inext:ifin] = autladat[1].data['ein_a_ref'][copyme]

  ladat = ladat[:ifin]

  ladat2['data']=ladat
  autlafilename_out = re.sub('autos.fits','acx.fits',autlafilename)
  autlvfilename_out = re.sub('autos.fits','acx.fits',autlvfilename)

  pyatomdb.util.write_la_file(autlafilename_out, ladat2, clobber=True)
  print("file written: %s"%(autlafilename_out))

  rtrck = open(rtrckfname, 'a')
  rtrck.write('LA File written\n')
  rtrck.close()

  ilvlist = numpy.where(lvdat2['data']['arad_ref']=='FIXME')[0]
  for ilv in ilvlist:
    lvdat2['data']['arad_tot'][ilv] = numpy.sum( ladat2['data']['einstein_a'][ ladat2['data']['upper_lev']==ilv+1])
    lvdat2['data']['arad_ref'][ilv] = 'AUTOSTRUCTRE_CX'
    lvdat2['data']['aaut_ref'][ilv] = 'AUTOSTRUCTRE_CX'

  pyatomdb.util.write_lv_file(autlvfilename_out, lvdat2, clobber=True)
  print("file written: %s"%(autlvfilename_out))
  rtrck = open(rtrckfname, 'a')
  rtrck.write('LV File written\n')
  rtrck.close()

  return autlvfilename_out, autlafilename_out

def find_valence_shells(lvdat, filemap):
  """
  Find the since excitation levels from an ion

  PARAMETERS
  ----------
  lvdat : hdulist
    The energy level file from pyatomdb.atomdb.get_data

  RETURNS
  -------
  isvalence : array(bool)
    True if level is a single valence shell excitation, false otherwise.
  ncapture : array(int)
    n shell captured in to (0 if n/a)
  isvalence : array(bool)
    l shell captured in to (-1 if n/a)
  """

  Z = lvdat[1].header['ELEMENT']
  z1 = lvdat[1].header['ION_STAT']+1
  z1p = z1+1

  isvalence = numpy.zeros(len(lvdat[1].data), dtype=bool)
  ncapture = numpy.zeros(len(lvdat[1].data), dtype=int)
  lcapture = numpy.zeros(len(lvdat[1].data), dtype=int)
  lcapture[:] = -1



  shlmax = sum(range(1,max(lvdat[1].data['N_QUAN'])+1))

   # get the ground state
  if z1==Z:
    # hydrogenic - everything is 1 electron excitation by definition...
#    isvalence[:] = True
    groundocc = numpy.zeros(shlmax, dtype=int)
  else:
    plvdat = pyatomdb.atomdb.get_data(Z, z1p, 'LV', settings = {'filemap': filemap})

    gndconfig = plvdat[1].data['ELEC_CONFIG'][0]
    pnel = Z+1-z1p

    groundocc, zzz = pyatomdb.atomic.config_to_occup(gndconfig, pnel, shlmax=shlmax)


  nlist = numpy.zeros(shlmax, dtype=int)
  llist = numpy.zeros(shlmax, dtype=int)
  nlist[0] = 1

  for i in range(1, shlmax):
    if llist[i-1] == nlist[i-1]-1:
      nlist[i] = nlist[i-1]+1
      llist[i] = 0
    else:
      nlist[i] = nlist[i-1]
      llist[i] = llist[i-1]+1

  nel = Z+1-z1
  # now go through each level
  for i in range(len(lvdat[1].data)):
    lvocc, zzz = pyatomdb.atomic.config_to_occup(lvdat[1].data['ELEC_CONFIG'][i], nel, shlmax=shlmax)
    if (sum(lvocc<groundocc)==0):
      if (sum(lvocc > groundocc)==1):
        # it's a valence shell!
        isvalence[i] = True
        ii = numpy.where(lvocc>groundocc)[0][0]
        ncapture[i] = nlist[ii]
        lcapture[i] = llist[ii]



  return isvalence, ncapture, lcapture



def make_cx_spectrumfile(Z, z1, donorsym, ftype, sigma_cx, sigma_cx_fname, fmapfile=False):
  """
  Create the necessary spectrum files for this ion for CX

  INPUTS
  ------
  Z : int
    atomic numbers
  z1 : int
    charge +1 of recombining ion
  donorsym : string
    Donor symbol
  ftype : string
    data type (nl, n, or none)
  sigma_cx : numpy.array
    The charge exchange data from calc_cross_sections
  sigma_cx_fname : string
    File name to write out to
  fmapfile : string
    Name of filemap file (STRONGLY RECOMMENDED: copy original, could
    be updated during processing)
  RETURNS
  -------
  filename : str
    The file where the spectrum was produced
  """

  z1r = z1-1 # charge of recombined ion
  settings = {'filemap':fmapfile,\
              'WriteIon':True,\
              'LinearGrid':True, \
              'GridMinimum':0.01, \
              'GridMaximum':100.0,\
              'NumGrid':100000,\
              'TwoPhoton':True,\
              'MinEpsilon':1e-10}
  linedata = numpy.zeros(0, dtype=generate_datatypes('linetype_cx_nl'))
  print("Hello, file type is %s"%(ftype))
  if ((ftype=='nl') | (ftype=='n')):
    # get the range of n or nl in the files

    if ftype=='n':

      nlist = pyatomdb.util.unique(sigma_cx['n'])
      nlist.sort()
      # get the capture into relevant l shell

      # ATOMIC DATA CHECK - does current data go to high enough n,l

      check_atomic_data(Z, z1, max(nlist), fmapfile=settings['filemap'])

      # get level populating mechanism

      # methodology: create a spectrum for each n shell & l distribution
      # captured in to (so n_shells * 4 l distributions)
      datacache={}
      lvdat = pyatomdb.atomdb.get_data(Z, z1r, 'LV', datacache=datacache, settings = settings)
      ladat = pyatomdb.atomdb.get_data(Z, z1r, 'LA', datacache=datacache, settings = settings)

      # get the list of n shells we need
      nmin = min(sigma_cx['n'])
      nmax = max(sigma_cx['n'])

      # find the levels which are valence (not inner shell)
      isvalence, nshell, lshell = find_valence_shells(lvdat, fmapfile)

      # assemble the pop matrix

      matrixA = numpy.zeros([len(lvdat[1].data), len(lvdat[1].data)],dtype=float)

      # element abundance
      abundances = pyatomdb.atomdb.get_abundance()
      abund=abundances[Z]

      ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', settings=settings, datacache=datacache)

      matrixA_in = {}
      matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
      pyatomdb.apec.gather_rates(Z, z1r, 1.0, 1.0, datacache=datacache, settings = settings,\
                   do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                   do_ir=False)

      for i in range(len(matrixA_in['init'])):
        matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]
      matrixA_in = {}

      for n in range(min(nlist),max(nlist)+1):
        for ldistname in ["even", "statistical", "landauzener","separable"]:
          ldist = calc_l_dstn(ldistname, n, z1r)

          capturerate= numpy.zeros(len(lvdat[1].data), dtype=float)
          continuum={}
          for l in range(len(ldist)):
            i = numpy.where((lshell == l) & (nshell == n))[0]
            if len(i) > 0:
              capturerate[i] = ldist[l]*lvdat[1].data['LEV_DEG'][i]/sum(lvdat[1].data['LEV_DEG'][i])
          print('capturerate:', capturerate)
          if sum(capturerate) > 0:
            capturerate*=-1

            levpop = pyatomdb.apec.calc_cascade_population(matrixA, capturerate)*abund

            linelist, tmptwophot = \
               pyatomdb.apec.do_lines(Z, z1r, levpop , 1.0, datacache=datacache, settings=settings, z1_drv_in=z1)
            continuum['twophot']= tmptwophot


            MinEpsilon = settings['MinEpsilon']
            ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                                      settings['GridMinimum'], \
                                      settings['GridMaximum'], \
                                      settings['NumGrid'])
            pseudocont = numpy.zeros(len(ebins)-1, dtype=float)
            print("nlines: %i, max %e"%(len(linelist), max(linelist['epsilon'])))
            if len(linelist) > 0:
              linelist=linelist[(linelist['lambda']>pyatomdb.const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                         (linelist['lambda']<pyatomdb.const.HC_IN_KEV_A /settings['GridMinimum'])]


            if len(linelist) > 0:
              weaklines = linelist[(linelist['epsilon']< MinEpsilon)]
              pseudocont,zzz = numpy.histogram(pyatomdb.const.HC_IN_KEV_A/weaklines['lambda'], bins=ebins, weights=weaklines['epsilon'])

              linelist = linelist[linelist['epsilon'] >= MinEpsilon]
            if settings['WriteIon']==True:
              ret = {}
              ret['cxdatatype']='kronos_n'
              ret['lines'] = linelist
              ret['continuum'] = continuum
              ret['pseudocont'] = pseudocont
              ret['settings'] = settings
              ret['abund'] = abund
              ret['n'] = n
              ret['ldistname'] = ldistname
              ret['Z'] = Z
              ret['z1'] = z1
              ret['donor'] = donorsym


              fname = "pickles/%s_%i_%i_n%i_l%s.pkl"%(donorsym, Z,z1,n,ldistname)
              pickle.dump(ret, open(fname, 'wb'))
              pickle_to_fits(fname, 'apec_cx2_%s_%i_%i_%i_%s'%(donorsym,Z, z1, n, ldistname))


    elif ftype=='nl':
      nlist = pyatomdb.util.unique(sigma_cx['n'])
      Slist = pyatomdb.util.unique(sigma_cx['s'])
      if len(Slist) > 1:
        Sresolved = True
      else:
        Sresolved = False


      check_atomic_data(Z, z1, max(nlist), fmapfile=fmapfile)

      # create the capture into each shell
      # methodology: create a spectrum for each n shell & l distribution
      # captured in to (so n_shells * 4 l distributions)
      datacache={}
      lvdat = pyatomdb.atomdb.get_data(Z, z1r, 'LV', datacache=datacache, settings = settings)
      ladat = pyatomdb.atomdb.get_data(Z, z1r, 'LA', datacache=datacache, settings = settings)

      # get the list of n shells we need
      nmin = min(sigma_cx['n'])
      nmax = max(sigma_cx['n'])

      # find the levels which are valence (not inner shell)
      isvalence, nshell, lshell = find_valence_shells(lvdat, fmapfile)

      # assemble the pop matrix

      matrixA = numpy.zeros([len(lvdat[1].data), len(lvdat[1].data)],dtype=float)

      # element abundance
      abundances = pyatomdb.atomdb.get_abundance()
      abund=abundances[Z]

      ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', settings=settings, datacache=datacache)

      matrixA_in = {}
      matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
      pyatomdb.apec.gather_rates(Z, z1r, 1.0, 1.0, datacache=datacache, settings = settings,\
                   do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                   do_ir=False)

      for i in range(len(matrixA_in['init'])):
        matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]
      matrixA_in = {}



##
      for n in range(min(nlist),max(nlist)+1):

        for l in range(0,n):
          for S in Slist:
            ldist = numpy.zeros(n, dtype=float)
            ldist[l] = 1.0

            capturerate= numpy.zeros(len(lvdat[1].data), dtype=float)
            continuum={}
            Stmp = numpy.array(lvdat[1].data['S_QUAN']*2+1, dtype=int)
            Stest = int(S*2+1)
            i = numpy.where((lshell == l) & (nshell == n) &(Stmp==Stest))[0]
            if len(i) > 0:
              capturerate[i] = ldist[l]*lvdat[1].data['LEV_DEG'][i]/sum(lvdat[1].data['LEV_DEG'][i])
            if sum(capturerate) > 0:
              capturerate*=-1*abund

              levpop = pyatomdb.apec.calc_cascade_population(matrixA, capturerate)


              linelist, tmptwophot = \
                 pyatomdb.apec.do_lines(Z, z1r, levpop , 1.0, datacache=datacache, settings=settings, z1_drv_in=z1)
              continuum['twophot']= tmptwophot

##

              MinEpsilon = settings['MinEpsilon']
              ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                                        settings['GridMinimum'], \
                                        settings['GridMaximum'], \
                                        settings['NumGrid'])
              pseudocont = numpy.zeros(len(ebins)-1, dtype=float)





              if len(linelist) > 0:
                linelist = linelist[(linelist['lambda']>pyatomdb.const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                                    (linelist['lambda']<pyatomdb.const.HC_IN_KEV_A /settings['GridMinimum'])]

              if len(linelist) > 0:
                weaklines = linelist[(linelist['epsilon']< MinEpsilon)]

                print("filtered out % i weak lines, left %i strong lines"%(len(weaklines), len(linelist)))

                pseudocont,zzz = numpy.histogram(pyatomdb.const.HC_IN_KEV_A/weaklines['lambda'], bins=ebins, weights=weaklines['epsilon'])
#      for iline, line in enumerate(weaklines):
#        e = pyatomdb.const.HC_IN_KEV_A /line['lambda']
#        ibin = numpy.where(ebins>e)[0][0] - 1
#        if (iline+1)%100000==0:
#          print "  on line %i"%(iline)
#        pseudocont[ibin]+=line['epsilon']

                linelist = linelist[linelist['epsilon'] >= MinEpsilon]
              print("Finish filtering linelist at %s"%(time.asctime()))


              print("post-filter: nlines: %i"%(len(linelist)))
              if settings['WriteIon']==True:
                ret = {}
                ret['cxdatatype']='kronos_nl'
                ret['lines'] = linelist
                ret['continuum'] = continuum
                ret['pseudocont'] = pseudocont
                ret['settings'] = settings
                ret['abund'] = abund
                ret['n'] = n
                ret['l'] = l
                ret['S'] = S
                ret['Z'] = Z
                ret['z1'] = z1
                ret['donor'] = donorsym


                fname = "pickles/%s_%i_%i_n%i_l%i_S%.1f.pkl"%(donorsym, Z,z1,n,l,S)
                pickle.dump(ret, open(fname, 'wb'))
                pickle_to_fits(fname, 'apec_cx2_%s_%i_%i_%i_%i_%.1f'%(donorsym, Z, z1, n, l, S))

#      for iE in range(sigma_cx['nC'][0]):
#        capturerate= numpy.zeros(len(lvdat[1].data), dtype=float)
#        continuum={}
#        print(iE,sigma_cx['nC'][0] )
#        for s in sigma_cx:
#          if s['s'] >= 0:
#            i = numpy.where((lshell == s['l']) & (nshell == s['n']) & (isvalence==True) & (lvdat[1].data['s_quan']==s['s']))[0]
#          else:
#            i = numpy.where((lshell == s['l']) & (nshell == s['n']) & (isvalence==True))[0]

#          if len(i) > 0:
#            capturerate[i] = s['C'][iE] * lvdat[1].data['LEV_DEG'][i]/sum(lvdat[1].data['LEV_DEG'][i])

#        if sum(capturerate) > 0:
#          capturerate*=-1*abund

#          levpop = pyatomdb.apec.calc_cascade_population(matrixA, capturerate)


#          linelist, tmptwophot = \
#            pyatomdb.apec.do_lines(Z, z1r, levpop , 1.0, datacache=datacache, settings=settings, z1_drv_in=z1)
#          continuum['twophot']= tmptwophot


          # MinEpsilon = settings['MinEpsilon']
          # ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                                    # settings['GridMinimum'], \
                                    # settings['GridMaximum'], \
                                    # settings['NumGrid'])
          # pseudocont = numpy.zeros(len(ebins)-1, dtype=float)
          # print("nlines: %i, max %e"%(len(linelist), max(linelist['epsilon'])))
          # if len(linelist) > 0:
            # weaklines = linelist[(linelist['epsilon']< MinEpsilon) &\
                       # (linelist['lambda']>pyatomdb.const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                       # (linelist['lambda']<pyatomdb.const.HC_IN_KEV_A /settings['GridMinimum'])]

            # for line in weaklines:
              # e = pyatomdb.const.HC_IN_KEV_A /line['lambda']
              # ibin = numpy.where(ebins>e)[0][0] - 1
              # pseudocont[ibin]+=line['epsilon']
            # linelist = linelist[linelist['epsilon'] >= MinEpsilon]
          # if settings['WriteIon']==True:
            # ret = {}
            # ret['lines'] = linelist
            # ret['continuum'] = continuum
            # ret['pseudocont'] = pseudocont
            # ret['settings'] = settings
            # ret['abund'] = abund
            # ret['iE'] =[iE]
            # print(sigma_cx['E'].shape)
            # ret['energy'] = sigma_cx['E'][0,iE]


            # fname = "pickles/%i_%i_iE%i.pkl"%(Z,z1,iE)
            # pickle.dump(ret, open(fname, 'wb'))

  elif ftype=='None':
    # In here, obtain maximum n from the formulae of doom.
    q = z1*1.0 # charge of recombin*ing* ion, == z1 of recombined ion
    Ip_d_list = {}
    Ip_d_list['h'] = 13.59843449
    Ip_d_list['he'] = 24.58738880

    if donorsym.lower() in Ip_d_list.keys():
      Ip_d = Ip_d_list[donorsym.lower()]
    else:
      print("Unknown donor %s: please add to make_cx_spectrum"%(donorsym))
      exit()

    nmax = int(numpy.ceil(q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                          (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))))

    print("nmax=", nmax)
    check_atomic_data(Z, z1, nmax, fmapfile=fmapfile)



    # get level populating mechanism

    # methodology: create a spectrum for each n shell & l distribution
    # captured in to (so n_shells * 4 l distributions)
    datacache={}
    lvdat = pyatomdb.atomdb.get_data(Z, z1r, 'LV', datacache=datacache, settings = settings)
#    ladat = pyatomdb.atomdb.get_data(Z, z1r, 'LA', datacache=datacache, settings = settings)


    # find the levels which are valence (not inner shell)
    isvalence, nshell, lshell = find_valence_shells(lvdat, fmapfile)

    # assemble the pop matrix
    matrixA = numpy.zeros([len(lvdat[1].data), len(lvdat[1].data)],dtype=float)

    # element abundance
    abundances = pyatomdb.atomdb.get_abundance()
    abund=abundances[Z]

    ladat = pyatomdb.atomdb.get_data(Z, z1, 'LA', settings=settings, datacache=datacache)

    matrixA_in = {}
    matrixA_in['init'], matrixA_in['final'], matrixA_in['rate']=\
    pyatomdb.apec.gather_rates(Z, z1r, 1.0, 1.0, datacache=datacache, settings = settings,\
                 do_la=True, do_ai=False, do_ec=False, do_pc=False,\
                 do_ir=False)

    for i in range(len(matrixA_in['init'])):
      matrixA[matrixA_in['final'][i], matrixA_in['init'][i]]+=matrixA_in['rate'][i]
    matrixA_in = {}

#HEREHEREHERE

    nprime = q *( ((1000*pyatomdb.const.RYDBERG/Ip_d)**0.5) * \
                        (( 1 + (q-1.0)/(2*q)**0.5)**-0.5))


    for n in [int(numpy.floor(nprime)), int(numpy.ceil(nprime))]:
      print("DOING n=%i, nprime=%f"%(n, nprime))

      for ldistname in ["even", "statistical", "landauzener","separable"]:
        print("starting ldistname = %s"%(ldistname))
        ldist = calc_l_dstn(ldistname, n, z1r)
        print("ldist:", ldist)
        capturerate= numpy.zeros(len(lvdat[1].data), dtype=float)
        continuum={}
        for l in range(len(ldist)):
          i = numpy.where((lshell == l) & (nshell == n))[0]
          print("There are %i levels with n=%i, l=%i. ldist[%i]=%e, sumLEVDEG=%i"%\
                (len(i),n,l,l,ldist[l], sum(lvdat[1].data['LEV_DEG'][i])))
          if len(i) > 0:
            capturerate[i] = ldist[l]*lvdat[1].data['LEV_DEG'][i]*1.0/sum(lvdat[1].data['LEV_DEG'][i])
        print("Capturerate = %e"%(sum(capturerate)), capturerate)
        if sum(capturerate) > 0:
          capturerate*=-1*abund

          levpop = pyatomdb.apec.calc_cascade_population(matrixA, capturerate)

          print("calling apec_do_lines: Z=%i, z1r=%i, z1_drv_in=%i"%(Z, z1r, z1))
          linelist, tmptwophot = \
             pyatomdb.apec.do_lines(Z, z1r, levpop , 1.0, datacache=datacache, settings=settings, z1_drv_in=z1)
          continuum['twophot']= tmptwophot


          MinEpsilon = settings['MinEpsilon']
          ebins = pyatomdb.apec.make_vector_nbins(settings['LinearGrid'], \
                                    settings['GridMinimum'], \
                                    settings['GridMaximum'], \
                                    settings['NumGrid'])
          pseudocont = numpy.zeros(len(ebins)-1, dtype=float)
          print("nlines: %i, max %e"%(len(linelist), max(linelist['epsilon'])))

          if len(linelist) > 0:

            linelist=linelist[(linelist['lambda']>pyatomdb.const.HC_IN_KEV_A /settings['GridMaximum'])   &\
                       (linelist['lambda']<pyatomdb.const.HC_IN_KEV_A /settings['GridMinimum'])]


          if len(linelist) > 0:
            weaklines = linelist[(linelist['epsilon']< MinEpsilon)]
            pseudocont,zzz = numpy.histogram(pyatomdb.const.HC_IN_KEV_A/weaklines['lambda'], bins=ebins, weights=weaklines['epsilon'])

            linelist = linelist[linelist['epsilon'] >= MinEpsilon]


          if settings['WriteIon']==True:
            ret = {}
            ret['cxdatatype']='acx'

            ret['lines'] = linelist
            ret['continuum'] = continuum
            ret['pseudocont'] = pseudocont
            ret['settings'] = settings
            ret['abund'] = abund
            ret['n'] = n
            ret['ldistname'] = ldistname
            ret['Z'] = Z
            ret['z1'] = z1
            ret['donor'] = donorsym

            fname = "pickles/%s_%i_%i_n%i_l%s_acx1.pkl"%(donorsym, Z,z1,n,ldistname)
            pickle.dump(ret, open(fname, 'wb'))

            pickle_to_fits(fname, 'apec_cx2_%s_%i_%i_%i_%s'%(donorsym,Z, z1, n, ldistname))






def pickle_to_fits(picklefilename, fitsfilestem):
  """
  Assemble the fits files from the pickle files.

  PARAMETERS
  ----------

  picklefilename : string
    The name of the pickle file generated by CX to read
  fitsfilestem : string
    The stem of the fits file. Existing files will be appended, not overwritten.
    formats are fitsfilestem_cxline.fits and fitsfilestem_cxcoco.fits

  RETURNS
  -------
  None
  """

  # check if pickle file exists

  try:
    cxdat = numpy.load(picklefilename,allow_pickle=True)

  except:
    print("Supplied pickle file does not exist: %s"%(picklefilename))
    raise


  # test if fits files exist
  linefilename = fitsfilestem+"_cxline.fits"
  contfilename = fitsfilestem+"_cxcont.fits"

  if os.path.isfile(linefilename) & os.path.isfile(contfilename):
    fits_cx_line = pyatomdb.pyfits.open(linefilename, mode='update')
    fits_cx_cont = pyatomdb.pyfits.open(contfilename, mode='update')

    datamode = 'update'
  else:
    datamode = 'new'


  # generate the HDU for this data

  Z=cxdat['Z']
  z1=cxdat['z1']

  ebins = pyatomdb.apec.make_vector_nbins(cxdat['settings']['LinearGrid'], \
                                          cxdat['settings']['GridMinimum'], \
                                          cxdat['settings']['GridMaximum'], \
                                          cxdat['settings']['NumGrid'])



  econt, contin = pyatomdb.apec.compress_continuum(ebins, cxdat['continuum']['twophot'], pyatomdb.const.TOLERANCE, minval=1e-38)
  cxdat['pseudocont'][cxdat['pseudocont']<0] =0
  epseudo, pseudocont = pyatomdb.apec.compress_continuum(ebins, cxdat['pseudocont'], pyatomdb.const.TOLERANCE, minval=1e-38)

  dat = {}
  dat['cont'] = numpy.zeros(1, dtype=pyatomdb.apec.generate_datatypes('continuum', npseudo=len(pseudocont), ncontinuum=len(contin)))

  dat['cont']['Z'][0] = Z
  dat['cont']['rmJ'][0] = z1
  dat['cont']['N_Cont'][0] = len(contin)
  dat['cont']['E_Cont'][0][:len(contin)] = econt
  dat['cont']['Continuum'][0][:len(contin)] = contin

  dat['cont']['N_Pseudo'][0] = len(pseudocont)
#  print('peudocont', pseudocont)
#  print('E_peudocont', epseudo)

  dat['cont']['E_Pseudo'][0][:len(pseudocont)] = epseudo
  dat['cont']['Pseudo'][0][:len(pseudocont)] = pseudocont


  linedata = numpy.zeros(0,dtype=generate_datatypes('linetype'))
  cocodata = numpy.zeros(0,dtype=generate_datatypes('continuum', ncontinuum=0, npseudo=0))

  linedata = numpy.append(linedata, cxdat['lines'])
  cocodata = pyatomdb.apec.continuum_append(cocodata, dat['cont'])



      # now make an HDU for all of this
  print(linedata.dtype.names)
  LHDUdat = pyatomdb.apec.create_lhdu_nei(linedata)


#      iseclHDUdat=iDens+iTe*settings['NumDens']
  LHDUdat.header['EXTNAME']=("CX_EMISS","name of this binary table extension")
#  LHDUdat.header['EXTVER']=(iseclHDUdat+1,"Index for this EMISSIVITY extension")

  LHDUdat.header['HDUNAME'] = ("CX DATA",\
                             'Spectral emission data')
  LHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
  LHDUdat.header['HDUCLAS1']=("CXLINE MODEL",\
                             'Line emission spectral model')
  LHDUdat.header['HDUCLAS2']=("CXLINE",\
                             'Charge Exchange Emission line data')

  LHDUdat.header['HDUVERS1']=("2.0.0",\
                               'version of format')

  LHDUdat.header['Z']=(Z,\
                             'Atomic number of recombining ion')

  LHDUdat.header['Z1']=(z1,\
                             'Ion charge +1 of recombining ion')

  LHDUdat.header['DONOR']=(cxdat['donor'],\
                             'Donor')


  LHDUdat.header['CXTYPE']=(cxdat['cxdatatype'],\
                             'Charge Exchange Data Type')
  LHDUdat.header['N']=(cxdat['n'],\
                             'N shell captured into')

  if cxdat['cxdatatype'] in ['kronos_n', 'acx']:
    LHDUdat.header['L']=(cxdat['ldistname'],\
                               'l shell distribution')
  else:
    LHDUdat.header['L']=(cxdat['l'],\
                               'l shell captured into')
    LHDUdat.header['S']=(cxdat['l'],\
                               'S Spin captured into')



  lhdulist= [LHDUdat]

  PrilHDU = pyatomdb.pyfits.PrimaryHDU()
  lhdulist.insert(0,PrilHDU)

  tmplhdulist = pyatomdb.pyfits.HDUList(lhdulist)
  tmplhdulist.writeto('fitstmp/%s'%(linefilename), overwrite=True, checksum=True)




  #now repeat for continuum data



      # now make an HDU for all of this
  CHDUdat = pyatomdb.apec.create_chdu_cie(cocodata)


#      isecCHDUdat=iDens+iTe*settings['NumDens']
  CHDUdat.header['EXTNAME']=("CX_EMISS","name of this binary table extension")
#  CHDUdat.header['EXTVER']=(isecCHDUdat+1,"Index for this EMISSIVITY extension")

  CHDUdat.header['HDUNAME'] = ("CX DATA",\
                             'Spectral emission data')
  CHDUdat.header['HDUCLASS'] = ("Proposed OGIP",\
                             'Proposed OGIP standard')
  CHDUdat.header['HDUCLAS1']=("CXLINE MODEL",\
                             'Continuum emission spectral model')
  CHDUdat.header['HDUCLAS2']=("CXCONT",\
                             'Charge Exchange Emission continuum data')

  CHDUdat.header['HDUVERS1']=("2.0.0",\
                               'version of format')

  CHDUdat.header['Z']=(Z,\
                             'Atomic number of recombining ion')

  CHDUdat.header['Z1']=(z1,\
                             'Ion charge +1 of recombining ion')

  CHDUdat.header['DONOR']=(cxdat['donor'],\
                             'Donor')


  CHDUdat.header['CXTYPE']=(cxdat['cxdatatype'],\
                             'Charge Exchange Data Type')
  CHDUdat.header['N']=(cxdat['n'],\
                             'N shell captured into')

  if cxdat['cxdatatype'] in ['kronos_n', 'acx']:
    CHDUdat.header['L']=(cxdat['ldistname'],\
                               'l shell distribution')
  else:
    CHDUdat.header['L']=(cxdat['l'],\
                               'l shell captured into')
    CHDUdat.header['S']=(cxdat['l'],\
                               'S Spin captured into')


  CHDUlist= [CHDUdat]

  PriCHDU = pyatomdb.pyfits.PrimaryHDU()
  CHDUlist.insert(0,PriCHDU)

  tmpCHDUlist = pyatomdb.pyfits.HDUList(CHDUlist)
  tmpCHDUlist.writeto('fitstmp/%s'%(contfilename), overwrite=True, checksum=True)




  #if settings['Ionization']=='CIE':
  #  tmpchdulist.writeto('%s_coco.fits'%(fileroot), clobber=True, checksum=True)
  #elif settings['Ionization']=='NEI':
  #  tmpchdulist.writeto('%s_comp.fits'%(fileroot), clobber=True, checksum=True)





# START HERE
if __name__=='__main__':

  # FIND THE ION
  ### z1 is the recombining ion charge +1

  Z = int(sys.argv[1])
  z1 = int(sys.argv[2])
  donorsym = sys.argv[3]

  # FIND WHICH DATA EXISTS

  # check for Kronos n,l resolved data
  fname, ftype = read_kronos_filedata(Z, z1, donorsym)

  sigma_cx_fname = 'acxdatabase/%i_%i_%s_CX_v1_0_0_a.fits'%(Z, z1, donorsym)
  print("ftype=", ftype)
  if ftype == 'nl':
    # nl resolved data
    sigma_cx = calc_cross_sections(Z, z1, fname, ftype)
  elif ftype == 'n':
    # n resolved data
    sigma_cx = calc_cross_sections(Z, z1, fname, ftype)
  elif ftype == 'None':
    sigma_cx = False
    # no data
    pass
  else:
    pass

  # store the cx cross sections
  if not ftype=='None':
    write_sigma_cx(sigma_cx_fname, sigma_cx, Z, z1, donorsym, ftype, fname)



  # Now we need to make a spectrum file
  ###HYDRACHANGE
  make_cx_spectrumfile(Z, z1, donorsym, ftype, sigma_cx, sigma_cx_fname, fmapfile='/export1/projects/ACX2/generation/filemap_3.0.9_cx')
#  make_cx_spectrumfile(Z, z1, donorsym, ftype, sigma_cx, sigma_cx_fname, fmapfile='/export1/projects/atomdb_cx/filemap_v3.0.9_cx')

  #

#if __name__=='__main__':
  #Z = int(sys.argv[1])
  #z1 = int(sys.argv[2])
##  acxtype = int(sys.argv[3])
  #donor = sys.argv[3]
  #datacache={}
  #for acxtype in range(1,17):
    #fname = 'results/acx_%i_donor_%s_Z_%i_rmJ_%i.pkl'%(acxtype,donor,Z,z1)
    #if os.path.exists(fname):
      #try:
        #t = numpy.load(fname)
      #except IOError:
        ##t = open(fname, 'r')
        ##tt=t.readline()
        ##t.close()
        ##if tt=='started\n':
        ##  os.remove(fname)
        #pass
    #if not os.path.exists(fname):
      #o = open(fname,'w')
      #o.write('started\n')
      #o.close()
      #sigma = process_ion(Z,z1, acxtype, donor, datacache=datacache)

      #if sum(sigma[1:]) > 0:
        #spectrum = create_spectrum(Z, z1,acxtype, donor, sigma, datacache=datacache)
      #else:
        #o = open(fname,'w')
        #o.write('DeliberatelyIgnore\n')
        #o.close()

  print("EXITING CLEANLY")

