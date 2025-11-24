from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(["select_recoil.py", "scatter.py", "trajectory.py", "estop.py"], language_level=3)
)