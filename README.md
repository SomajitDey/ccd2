[![CI Build Status](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml/badge.svg)](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml)

# Collective Cell Dynamics
This software contains the simulation engine and analysis tools for the bead-spring model as presented in <u>A. Mkrtchyan, J. Astrom, M. Karttunen, *Soft Matter*, 2014, **10**, 4332</u>. The model may or may not have been modified to suit our needs.

# Status
Work in progress

# Dependencies
- [fortdepend](https://github.com/ZedThree/fort_depend.py) or its [fork](https://github.com/PhyBi/fortdepend) for generating dependencies during build
- [pv](https://command-not-found.com/pv) for showing real-time progress bar

# Build
Use `make` as usual. Other uses: `make rebuild`, `make clean`.
Default compiler: `gfortran`. To use `ifort` instead, `make FC=ifort`
To include OpenMP, set the `OMP` variable when running `make`, e.g. `make OMP=set rebuild`.

# Install
`make install`

If you don't want to install it right away but test: `. setup_test_env.sh`

# Uninstall
`make uninstall`

# Command-line
`ccd`

# Sample Workflow
Edit the `params.in` file. Put only those key-value pairs that you don't want the default for.

```bash
# Random initialization
ccd rinit

# Run
ccd run [-a | --append] [-n | --no-status-dump] [-f | --force]

# Checkpoint to XY dump
ccd cpt_to_xy

# Visualize using gnuplot
ccd visual

# To check live run progress
ccd status
```

# Note
If met with segmentation faults or stack-smashing error, make the stack size unlimited in the Bash session with `ulimit -s unlimited`.
If the problem persists, rebuild with `make DEBUG=set` and report @ [issue](https://github.com/PhyBi/Collective-Cell-Dynamics/issues).
Also try the `-heap-arrays` flag (FF for ifort) in the Makefile.

# License
No use of this software shall be made without permission from the PI, [Dr. Dipjyoti Das](mailto:dipjyoti.das@iiserkol.ac.in).
