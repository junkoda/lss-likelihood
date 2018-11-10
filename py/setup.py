#
# $ make
#

from distutils.core import setup, Extension
#import numpy as np
import os


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# directories for include -I(idir)

idirs = []
if 'IDIRS' in os.environ:
    idirs = os.environ('IDIRS')
    if idirs:
        idirs = idirs.split()

#idirs = ['../lib', np.get_include()] #+ idirs


# directories for libraries -L(dir)

ldirs = []
if 'LDIRS' in os.environ:
    ldirs = os.environ["LDIRS"]
    if ldirs:
        ldirs = ldirs.split()



libs = ['gsl'] #os.environ['LIBS'].split()

setup(name='lss_likelihood',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['lss_likelihood.model',
                  'lss_likelihood.covariance',
                  'lss_likelihood.data',
      ],
      ext_modules=[
          Extension('lss_likelihood._lss_likelihood',
                    ['py_package.cpp', 'py_util.cpp',
                     'py_model.cpp', 'py_data.cpp',
                     'py_likelihood.cpp',
                    ],
                    include_dirs = idirs,
                    libraries = libs,
                    undef_macros = ['NDEBUG'],
                    #extra_compile_args = [os.environ['OPT']],
                    #library_dirs = ldirs,

          )
      ],
      packages=['lss_likelihood'],
)
