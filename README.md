# This is a fork of LAMMPS extended by AutoPas

## How to run with AutoPas

### Compile

```bash
mkdir lammps-autopas/build && cd "$_"

cmake \
    -DPKG_USER-AUTOPAS=yes \
    -DPKG_USER-OMP=yes \
    -DPKG_KOKKOS=yes \
    -DKOKKOS_ARCH=HSW \
    -DKOKKOS_ENABLE_OPENMP=yes \
    -DBUILD_OMP=yes \
    -DWITH_GZIP=false \
    -DWITH_JPEG=false \
    -DWITH_PNG=false \
    ../cmake

make
```

### Run

Examples:
```bash
lmp -i in.lj -autopas on log debug -sf autopas -v t 1000 -v x 128 -v y 128 -v z 128 

```
or

```bash
lmp -i in.lj -autopas on notune t balanced-sliced-verlet c VerletListsCells d AoS n disabled estimator none -sf autopas -v t 1000 -v x 128 -v y 128 -v z 128 
```

## Original LAMMPS README

This is the LAMMPS software package.

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel
Simulator.

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

LAMMPS is a classical molecular dynamics simulation code designed to
run efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary author of the code is Steve Plimpton, who can be emailed
at sjplimp@sandia.gov.  The LAMMPS WWW Site at lammps.sandia.gov has
more information about the code and its uses.

The LAMMPS distribution includes the following files and directories:

* README.md			   this file
* LICENSE			   the GNU General Public License (GPL)
* bench			   benchmark problems
* cmake			   CMake build system
* doc			   documentation
* examples		   simple test problems
* lib			   libraries LAMMPS can be linked with
* potentials		   interatomic potential files
* python			   Python wrapper on LAMMPS as a library
* src			   source files
* tools			   pre- and post-processing tools

Point your browser at any of these files to get started:

* http://lammps.sandia.gov/doc/Manual.html         the LAMMPS manual
* http://lammps.sandia.gov/doc/Intro.html          hi-level introduction
* http://lammps.sandia.gov/doc/Build.html          how to build LAMMPS
* http://lammps.sandia.gov/doc/Run_head.html       how to run LAMMPS
* http://lammps.sandia.gov/doc/Developer.pdf       LAMMPS developer guide

You can also create these doc pages locally:

```bash
cd doc
make html                # creates HTML pages in doc/html
make pdf                 # creates Manual.pdf and Developer.pdf
```
