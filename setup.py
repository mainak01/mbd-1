#!/usr/bin/env python
from setuptools import setup
import setuptools.command.build_ext
import sys
if sys.version_info[0] == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser

config = ConfigParser()
with open('setup.cfg') as f:
    config.readfp(f)


class build_ext(setuptools.command.build_ext.build_ext):
    def build_extensions(self):
        self.compiler.linker_so[0] = config.get('fortran', 'mpifc')
        setuptools.command.build_ext.build_ext.build_extensions(self)


setup(
    name='pymbd',
    version='0.2.dev1',
    url='https://github.com/azag0/mbd',
    author='Jan Hermann',
    author_email='jh@janhermann.cz',
    license='MIT License',
    py_modules=['pymbd'],
    cffi_modules=['mbd_build.py:ffi'],
    setup_requires=['cffi>=1.4.2'],
    install_requires=['cffi>=1.4.2'],
    cmdclass={'build_ext': build_ext}
)
