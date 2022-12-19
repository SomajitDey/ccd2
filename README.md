[![CI Build Status](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml/badge.svg)](https://github.com/PhyBi/Collective-Cell-Dynamics/actions/workflows/build.yml)

# Collective Cell Dynamics
This software contains the simulation engine and analysis tools for the bead-spring model as presented in <u>A. Mkrtchyan, J. Astrom, M. Karttunen, *Soft Matter*, 2014, **10**, 4332</u>. The model may or may not have been modified to suit our needs.

# Status
Work in progress. Check [TODO](#TODO) list below.

# Dependencies
- [fortdepend](https://github.com/ZedThree/fort_depend.py) or its [fork](https://github.com/PhyBi/fortdepend) for generating dependencies during build
- [bash](https://command-not-found.com/bash) because we suck at writing fully POSIX-compliant shell scripts
- [pv](https://command-not-found.com/pv) for showing real-time progress bar
- [gnuplot](https://command-not-found.com/gnuplot) for visualization

# Build
- Install, if non-existent, the above [dependencies](#dependencies) first

- Download this project: `git clone --depth=1 https://github.com/PhyBi/Collective-Cell-Dynamics ccd`

- Change to the downloaded project directory: `cd ccd`

- Compile: Use `make` as usual. Other uses: `make rebuild`, `make clean`.

Notes: Default compiler is `gfortran`. To use `ifort` instead, use `make FC=ifort`. To include OpenMP, set the `OMP` variable when running `make`, e.g. `make OMP=set rebuild`. Set the `DEBUG` variable in order to build in debug mode.

# Install
`make install`

If you don't want to install it right away but test: `. setup_test_env.sh`

# Uninstall
`make uninstall`

# Command-line
`ccd`

# Sample Workflow
- Create and change into a new directory.

- Copy and edit the [params.in](/params.in) file. *Put only those key-value pairs for the run parameters that you don't want the default for*. `ccd show_params` lists all parameters ccd takes. For meaning and default values of the parameters, please look up [src/mod_parameters.f90](/src/mod_parameters.f90).

```bash
# Random initialization
ccd init

# Run
ccd run [-a | --append] [-n | --no-status-dump] [-f | --force] > metadata.txt 2> logfile.txt

# Checkpoint to XY dump
ccd cpt_to_xy

# Visualize using gnuplot
ccd visual

# To check live run progress
ccd status

# To archive run results (trajectory: traj.bin and final-state/checkpoint: state.cpt)
ccd archive

# To extract archived run results
ccd archive metadata.txt
```

- If you want to use another path instead of `params.in` for the run-parameter file, put the same in the enviroment variable `CCD_PARAMS_PATH`

- If all (or most) of your runs use certain common *non-default* parameter values, provide those in the key-value file (RC file) at `${HOME}/.ccdrc`. `ccd` reads these parameters before reading parameters from `params.in` or `${CCD_PARAMS_PATH}`. Values read in from the latter file take precedence in case of conflicts. Non-default locations of the RC file may be passed through the `CCD_RC_PATH` enviroment variable.

- It is also possible to pass run-time parameters using the command-line. These values take precedence over those read in from from `params.in` or `${CCD_PARAMS_PATH}` in case of conflicts. For example,
```bash
ccd -p '<parameterA>=<valueA>' -p '<parameterB>=<valueB>' show_params
```

# Note
If met with segmentation faults or stack-smashing error, make the stack size unlimited in the Bash session with `ulimit -s unlimited; export OMP_STACKSIZE=500m`.
If the problem persists, rebuild with `make DEBUG=set` and report @ [issue](https://github.com/PhyBi/Collective-Cell-Dynamics/issues).
Also try the `-heap-arrays` compiler flag (`FF`) for ifort in the Makefile.

# License
No use of this software shall be made without permission from the PI, [Dr. Dipjyoti Das](mailto:dipjyoti.das@iiserkol.ac.in).

# Disclaimer
Use this software at your own risk. We the devs or this organization/lab/institute or the PI are/is not responsible for any damage caused by the use of this software.

# TODO
This software is built from a monolithic legacy code. Hence much had and still would have to be done towards enhancements in performance, maintainability, user-interface etc.
- [x] giving the entire codebase a permanent home online
- [x] fragmentation and build automation
- [x] replacing formatted trajectory dump with much faster and more compact unformatted and asynchronous I/O
- [x] include a parameters input file instead of using hardcoded values
- [x] replacing component based code with vector notation wherever possible for SIMD
- [x] remove all GOTO statments
- [x] making the run trajectory reproducible (at least for serial/single-threaded execution) by storing the PRNG seeds
- [x] include real-time progress bar
- [x] include performance calculation and dump
- [x] establishing separate run modes using command line options
- [x] replacing long hardcoded I/O paths with short filenames
- [x] inclusion of logging (at stderr)
- [x] checkpointing. Checkpoints serve dual purpose - run recovery and initialization
- [x] command-line autocompletion
- [x] very basic multithreading using OpenMP
- [x] continuous integration and git-hooks
- [ ] beads are stored in rows. Store them in columns instead for much performance improvement (as fortran is column major)
- [x] replace the overcomplicated neighborlist structure with a simple bead-based linked list for decreased overhead
- [x] replace the array-bound-based implementation of circular boundary conditions of beads within cells for better maintainability 
- [x] building the cell-cell neighborlist and dumping it in trajectory file in the most compressed way for later analysis such as hexatic order parameter
- [x] current initialization works only for the hardcoded system size. The fix (which would also make the system size assignable by the user) is ready for deployment but can only come after the neighborlist fix.
- [ ] linting
- [ ] performance oriented profiling and polishing
- [x] enhancing the driver code (`ccd`) as well as the bash-completion script
- [ ] include detailed [docs](docs/)
- [ ] include analysis tools (legacy or new)
