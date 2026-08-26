from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

'''  To use uspex_softmode_core module, run this setup first!
python setup_softmode.py build_ext --inplace
'''

extensions = [
    Extension(
        "uspex_softmode_core",
        ["uspex_softmode_core.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-ffast-math", "-fno-finite-math-only", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    ),
]

setup(
    name="uspex_softmode_core",
    ext_modules=cythonize(extensions, compiler_directives={
        'language_level': 3,
        'boundscheck': False,
        'wraparound': False,
        'cdivision': True,
    }),
)

extensions = [
    Extension(
        "uspex_fast_core",
        sources=["uspex_fast_core.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-march=native", "-fopenmp"],
        extra_link_args=["-fopenmp"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]

setup(
    name="uspex_fast",
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "nonecheck": False,
        },
    ),
)

extensions = [
    Extension(
        "uspex_rotation_core",
        ["uspex_rotation_core.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-ffast-math","-fno-finite-math-only", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    ),
]

setup(
    name="uspex_rotation_core",
    ext_modules=cythonize(extensions, compiler_directives={
        'language_level': 3,
        'boundscheck': False,
        'wraparound': False,
        'cdivision': True,
    }),
)
