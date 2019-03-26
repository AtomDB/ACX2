from setuptools import setup, Extension
import os



def get_version(relpath):
    """read version info from file without importing it"""
    from os.path import dirname, join
    for line in open(join(dirname(__file__), relpath)):
        if '__version__' in line:
            if '"' in line:
                return line.split('"')[1]
            elif "'" in line:
                return line.split("'")[1]


readme = open('README.txt', 'r')
README_TEXT = readme.read()
readme.close()



setup(name='acx2',
      version=get_version('__init__.py'),
      description='AtomDB Charge Exchange python library.',
      url='http://www.atomdb.org',
      author='Adam Foster',
      author_email='afoster@cfa.harvard.edu',
      license='Smithsonian',
      py_modules=['acx2'],
      classifiers=['Development Status :: 4 - Beta',\
                   'Environment :: Console',\
                   'Intended Audience :: Developers',\
                   'Intended Audience :: Education',\
                   'Intended Audience :: End Users/Desktop',\
                   'Intended Audience :: Science/Research',\
                   'Topic :: Scientific/Engineering :: Astronomy',\
                   'Topic :: Scientific/Engineering :: Physics',\
                   'Programming Language :: Python :: 3',\
                   'Operating System :: POSIX',\
		   'Operating System :: MacOS'],
      zip_safe=False,
      long_description = README_TEXT,\
      install_requires=[
      "pyatomdb"])

