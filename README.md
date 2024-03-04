[![CI Build Status](https://github.com/SomajitDey/ccd2/actions/workflows/build.yml/badge.svg)](https://github.com/SomajitDey/ccd2/actions/workflows/build.yml) ![Bash completion](https://img.shields.io/badge/Commandline%20Completion-Enabled-green) ![OS](https://img.shields.io/badge/Platform-Linux%2C%20MacOS(X)%2C%20WSL-blue)

# Collective Cell Dynamics 2.0
This repository contains major changes to our earlier [project](https://github.com/PhyBi/Collective-Cell-Dynamics).

# Status
Work in progress. Unstable and incomplete. Not fit for public use.

# Dependencies
- [gfortran](https://command-not-found.com/gfortran) or [ifort](https://gist.github.com/SomajitDey/aeb6eb4c8083185e06800e1ece4be1bd)(recommended). Note: ifort is free now, doesn't require a license anymore.
- [fortdepend](https://github.com/ZedThree/fort_depend.py) or its [fork](https://github.com/PhyBi/fortdepend) for generating dependencies during build. For Ubuntu 18.04, you may also have to install *importlib-metadata* with: `pip3 install importlib.metadata`
- [bash](https://command-not-found.com/bash) because we suck at writing fully POSIX-compliant shell scripts
- [pv](https://command-not-found.com/pv) for showing real-time progress bar
- [gnuplot](https://command-not-found.com/gnuplot) for visualization
- [xz](https://command-not-found.com/xz) for trajectory compression
- [ffmpeg](https://command-not-found.com/ffmpeg) for movie generation
- [jq](https://command-not-found.com/jq), [curl](https://command-not-found.com/curl) and [sponge](https://command-not-found.com/sponge) for quotes
- [helpdoc](https://github.com/somajitdey/helpdoc) for showing help/usage documentation
- [freud](https://freud.readthedocs.io/en/latest/gettingstarted/installation.html) for constructing periodic Voronoi tesselations. It may be installed with: `pip3 install freud-analysis`. (freud uses [Voro++](https://github.com/chr1shr/voro), a C++ library, with which we may replace freud later on).

# Build
- Install, if non-existent, the above [dependencies](#dependencies) first

- Download this project:
```bash
git clone --depth=1 https://github.com/SomajitDey/ccd2 ccd
```

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

Supported subcommands are available using bash-completion: `ccd <TAB><TAB>`

# Sample Workflow
The information in this section is anything but exhaustive; the subcommands are more versatile than apparent from what follows. Use the `-h` option in each subcommand to see its detailed usage. For any queries, feel free to ask [us](https://github.com/PhyBi/Collective-Cell-Dynamics/discussions).

- Create and change into a new directory.

- Copy and edit the [params.in](/params.in) file. *Put only those key-value pairs for the run parameters that you don't want the default for*. `ccd show_params` lists all parameters ccd takes. For meaning and default values of the parameters, please look up [src/mod_parameters.f90](/src/mod_parameters.f90).

- Simulations also produce/use/update two *non-portable* binary files: the trajectory (`traj.bin`) and a checkpoint (`state.cpt`). The checkpoint serves dual purpose: run recovery and initialization. It helps in run recovery by storing the last uncorrupted record/frame number in the trajectory file as well as the number of pending timesteps and the current timestep. Initialization is served by storing the complete memory representation of the coordinates at a timepoint (state).

```bash
# Random initialization
ccd init

# Run
ccd run > 'metadata.txt' 2> 'logfile.txt'

# State/Checkpoint(state.cpt) to XY dump(config.xy)
ccd cpt_to_xy

# Visualize config from XY dump (config.xy), highlighting the 10th cell in green
ccd visual -L -H '10@green#' # Use the GNUPLOT_PATH env variable if gnuplot is not in PATH

# To check live run progress
ccd status

# Make a movie (.mp4) of length=15 seconds out of the trajectory (traj.bin)
ccd traj_to_xy 'movie'
ccd movie -l 15 -H1# 'movie'

# To archive run results (trajectory: traj.bin and final-state/checkpoint: state.cpt)
ccd archive 'metadata.txt'

# To retrieve/restore archived run results
ccd retrieve 'metadata.txt'

# Garbage cleanup
ccd clean_archive 'metadata.txt'
```

### Other ways to pass parameters
- If you want to use another path instead of `params.in` for the run-parameter file, put the same in the enviroment variable `CCD_PARAMS_PATH`

- If all (or most) of your runs use certain common *non-default* parameter values, provide those in the key-value file (RC file) at `${HOME}/.ccdrc`. `ccd` reads these parameters before reading parameters from `params.in` or `${CCD_PARAMS_PATH}`. Values read in from the latter file take precedence in case of conflicts. Non-default locations of the RC file may be passed through the `CCD_RC_PATH` enviroment variable.

- It is also possible to pass run-time parameters using the command-line. These values take precedence over those read in from from `params.in` or `${CCD_PARAMS_PATH}` in case of conflicts. For example,
```bash
ccd -p '<parameterA>=<valueA>; <parameterB>=<valueB>' show_params
```
or equivalently,
```bash
ccd -p '<parameterA>=<valueA>' -p '<parameterB>=<valueB>' show_params
```

### Organization made easy
- Redirect the stdout of `ccd run` to metadata files as shown above.

- Archive the run results by providing the metadata file to `ccd archive`.

- Organize the small metadata files with commit messages and branches using `git`, and preserve them in cloud (e.g. GitHub).

- Retrieve the run results later, as necessary, by providing the corresponding metadata file to `ccd retrieve`.

- Backup the archive (`${HOME}/.ccd`) from time to time, copying the new files only.

### Quotes
To turn off the GROMACS-like quotes at the end of each command, set the `CCD_NO_QUOTES` enviroment variable. Do so if you mostly use `ccd` when offline.

### Extensibility
User can put her own routines (source code and scripts) in the [custom](/custom/) directory, provided the main sources and executable scripts are named as `ccd_<subcmd>.f90` and `ccd_<subcmd>` respectively. Such sources will be automatically built (installed) when one [builds](#build) (installs) the current project. Built executables and scripts can be invoked as `ccd <subcmd> [<options>]`.

# Note
If met with segmentation faults or stack-smashing error, make the stack size unlimited in the Bash session with `ulimit -s unlimited; export OMP_STACKSIZE=500m`.
If the problem persists, rebuild with `make DEBUG=set` and report @ [issue](https://github.com/PhyBi/Collective-Cell-Dynamics/issues).
Also try the `-heap-arrays` compiler flag (`FF`) for ifort in the Makefile.

# License
TBD. In its current state this software is not fit for public use.

# Disclaimer
Use this software at your own risk. The developer(s) is(are) not responsible for any damage caused by the use of this software.
