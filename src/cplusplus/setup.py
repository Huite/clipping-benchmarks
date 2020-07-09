from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


ext_modules = [
    Extension(
        "cpp_clipping",
        ["cpp_clipping.pyx", "clipping.cpp"],
        language='c++',
        extra_compile_args=["-O2"]
    )
]

setup(cmdclass = {'build_ext': build_ext}, ext_modules = ext_modules, include_dirs=[numpy.get_include()])
