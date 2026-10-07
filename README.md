# SFD (SFD_Calc_P4)
Code to compute Strain Functional Descriptors. This is approved for release under O5065. 
# Compilation and running
g++ compiler

make (creates Moment.5 binary)

make clean (removes .o files)

./Moment.5 (runs the binary)

Requires a params.txt file. The lines in the file are as follows:

1. Filename with atomic coordinates
2. Number of atoms
3. Cutoff for neighbors
4. Sigma of the Gaussian
5. Binary variable: 1 implies input is in the SPaSM binary format
6. Binary variable: 1 implies the weighting uses a reference configuration
7. Binary variable: 1 implied reading LAMMPS dump file (id type x y z)

# License
SFD is distributed as open source software available under a GPL3 license.


# Copyright 
© 2026. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
