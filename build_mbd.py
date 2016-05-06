import sys
import cffi
import json

with open('config.json') as f:
    env = json.load(f)['env']
header_file = sys.argv[1]
ffi = cffi.FFI()
ffi.set_source(
    '_mbd',
    '#include "{0}"'.format(header_file),
    extra_objects=[sys.argv[2]]
)
with open(header_file) as f:
    ffi.cdef(f.read())

if __name__ == '__main__':
    ffi.compile()
