#!/usr/bin/env python3
from pymbd import Context, Geometry


ctx = Context()
geom = Geometry.from_path('benzene.xyz')
crystal = Geometry.from_path('crystal.aims')
print(ctx)
print(geom)
print(crystal)
print(crystal.get_free_atoms())
ctx.do_rpa = True
