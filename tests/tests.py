#!/usr/bin/env python3
from pymbd import Context, Geometry, Damping


ctx = Context()
geom = Geometry.from_path('benzene.xyz')
crystal = Geometry.from_path('crystal.aims')
damping = Damping('fermi,dip', 'mbd@rsscs', 'pbe')
print(ctx)
print()
print(crystal)
print()
print(crystal.get_free_atoms())
print()
print(damping)
ctx.do_rpa = True
