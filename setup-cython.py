'''
Created on Feb 3, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# local imports
from setup import get_cython_extension_modules, setup_kwargs

def main():
    setup(ext_modules=get_cython_extension_modules(),
          cmdclass={'build_ext': build_ext},
          **setup_kwargs)

if __name__ == '__main__':
    main()