from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name="matrix_algos_cy",
    ext_modules=cythonize(
        "matrix_algos_cy.pyx",
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "nonecheck": False,
        },
    ),
    include_dirs=[np.get_include()],
)
