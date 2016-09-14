from distutils.core import setup, Extension
import numpy.distutils.misc_util
import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

module1 = Extension('FMM3D', sources = ['3DFMMmodule.cpp', '3DFMM.cpp']) #, extra_compile_args=['-fopenmp'], extra_link_args=['-lgomp']

setup (name = 'FMM3D',
        version = '1.0',
        description = 'Package for calculating seismic travel time over regular 3D grid using Fast Marching Method',
        ext_modules = [module1],
	   include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs())
