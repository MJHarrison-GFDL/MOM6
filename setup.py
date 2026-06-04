"""
"""
from numpy.distutils.core import setup,Extension

doclines = __doc__.split("\n")


regrid  = Extension(name = 'regrid',
                include_dirs = ['.mod',],
                library_dirs = ['.'],
                libraries = ['Regrid'],
                sources = ['regrid.f90','remap.f90'])


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyRegridding",
          version = '1.0.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Matthew Harrison",
          author_email = "matthew.harrison@noaa.gov",
          url = "none",
          license = 'CCL',
          platforms = ["any"],
          ext_modules = [regrid],
          )
