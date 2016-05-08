# `pymbd` — Many-body dispersion package

The many-body dispersion ([MBD](http://www.fhi-berlin.mpg.de/~tkatchen/MBD/)) method evaluates van der Waals interactions on top of density functional theory ([DFT](https://en.wikipedia.org/wiki/Density_functional_theory)).

This package provides an implementation of MBD as a Fortran 95 module with a C interface (Fortran 2003) and Python 2/3 bindings.

To build and test, run

```
./configure
make build test
```

Configuration options can be listed with

```
./configure --help
```
To install locally with `pip`, run

```
make install
```