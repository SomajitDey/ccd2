[![CI Build Status](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml/badge.svg)](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml)

# Collective Cell Dynamics
This software contains the simulation engine and analysis tools for the bead-spring model as presented in <u>A. Mkrtchyan, J. Astrom, M. Karttunen, *Soft Matter*, 2014, **10**, 4332</u>. The model may or may not have been modified to suit our needs.

# Status
Work in progress

# Dependencies
- [fortdepend](https://github.com/ZedThree/fort_depend.py) or its [fork](https://github.com/PhyBi/fortdepend)
- [pv](https://command-not-found.com/pv)

# Build
Use `make` as usual. Other uses: `make rebuild`, `make clean`, `make install`, `make uninstall`.
Default compiler: `gfortran`. To use `ifort` instead, `make FC=ifort`
To include OpenMP, set the `OMP` variable when running `make`, e.g. `make OMP=set rebuild`.

# License
No use of this software shall be made without permission from the PI, [Dr. Dipjyoti Das](mailto:dipjyoti.das@iiserkol.ac.in).
