# __all__=["acx2"]

import acx2
import pyatomdb
import curl
import os
import shutil
import time
import urllib
import wget

__version__="2.1.3"


#-------------------------------------------------------------------------------

def download_acx2_emissivity_files(adbroot, userid, version):

  """
  Download the AtomDB equilibrium emissivity files for AtomDB"

  This code will go to the AtomDB FTP site and download the necessary files.
  It will then unpack them into a directory adbroot. It will not
  overwrite existing files with the same md5sum (to avoid pointless updates)
  but it will not know this until it has downloaded and unzipped the main
  file.

  Parameters
  ----------

  adbroot : string
    The location to install the data. Typically should match $ATOMDB
  userid : string
    An 8 digit ID number. Usually passed as a string, but integer
    is also fine (provided it is all numbers)
  version : string
    The version string for the release, e.g. "3.0.2"

  Returns
  -------
  None
  """
  #
  # Version 0.1 Initail Release
  # Adam Foster 24th September 2015
  #
  import tarfile
  # set up remote file name
  fname = "acx2-v%s-data.tar.bz2"%(version)

  # set up temporary directory to hold data

  if adbroot[0] != '/':
    # is a relative path
    adbroot = "%s/%s"%(os.getcwd(), adbroot)

  pyatomdb.util.mkdir_p(adbroot)
  if adbroot[-1]=='/':
    tmpdir = adbroot+'installtmp'
  else:
    tmpdir = adbroot+'/installtmp'


  print("making directory %s"%(tmpdir))
  pyatomdb.util.mkdir_p(tmpdir)

  # get the files
  urllib.request.urlcleanup()
  fnameout = wget.download('%s/releases/%s'%(pyatomdb.const.FTPPATH,fname), out="%s/%s"%(tmpdir, fname))
  # collect user statistics if allowed.
  pyatomdb.util.record_upload(fname)

  #uncompress
  print("")
  print("Uncompressing ... ", end=' ')

  tf = tarfile.open(name=fnameout, mode='r:bz2')
  tf.extractall(path=tmpdir)
  print("Uncompressed")
  # copy the files
  dirname = 'acx2-v%s-data'%(version)
  for l in os.listdir('%s/%s'%(tmpdir, dirname)):
    print("moving %s/%s/%s to %s/%s"%(tmpdir, dirname, l, adbroot, l))
    shutil.move("%s/%s/%s"%(tmpdir, dirname, l), "%s/%s"%(adbroot, l))

  print("...done")

  shutil.rmtree(tmpdir)



# check initialization
on_rtd=os.environ.get('READTHEDOCS')=='True'
install='no'
if on_rtd:
  pass
else:
  if os.environ.get('ATOMDB')==None:
    print("ATOMDB environment variable is not set. Please set this to where your "+\
          "AtomDB files are stored. If you do not have these files, open python and "+\
          "call `import pyatomdb' to set up. Exiting.")
    exit()

    # issues here
  adbroot=os.environ.get('ATOMDB')
  READMEFILE = os.path.expandvars("$ATOMDB/README_ACX")
  if not os.path.exists(READMEFILE):
    print("ACX2 emissivity files not found. Installing now to $ATOMDB folder (%s)."%\
          os.path.expandvars('$ATOMDB'))
    install='yes'

  a=curl.Curl()
  version=a.get('%s/releases/LATEST_ACX'%(pyatomdb.const.FTPPATH))[:-1].decode(encoding='ascii')
  userprefs = pyatomdb.util.load_user_prefs(adbroot=adbroot)
  userid = userprefs['USERID']
  userprefs['LASTCXVERSIONCHECK'] = time.time()
  userprefs = pyatomdb.util.write_user_prefs(userprefs, adbroot=adbroot)

  if install=='yes':
      download_acx2_emissivity_files(os.environ.get('ATOMDB'), userid, version)

#  elif not os.path.exists(os.path.expandvars("$ATOMDB/userdata")):
#    install = ''
#    while not install in ['yes','no']:
#      install = util.question("ATOMDB Environment Variable Not Set. Set up AtomDB files now?", "yes",multichoice=["yes","no","info"])
#      if install=='info':
#        print("If yes, this will run the pyatomdb.util.initialize() script, which installs the necessary data files to run pyatomdb")

#    if install=='yes':
#      util.initialize()


  # now get the CX files
  #elif not os.path.exists(os.path.expandvars("$ATOMDB/userdata")):



