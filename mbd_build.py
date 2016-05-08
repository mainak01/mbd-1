import cffi
import sys
import os
if sys.version_info[0] == 2:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser


config = ConfigParser()
with open('setup.cfg') as f:
    config.readfp(f)
sourcefile = 'src/mbd.h'
ffi = cffi.FFI()
ffi.set_source(
    '_mbd_backend',
    '#include "{0}"'.format(os.path.abspath(sourcefile)),
    libraries=['mbd'],
    library_dirs=config.get('cffi', 'library-dirs').split(),
    extra_link_args=config.get('cffi', 'lapack').split()
)
with open(sourcefile) as f:
    ffi.cdef(f.read())
