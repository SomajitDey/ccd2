# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Fortran Compiler
FC := gfortran

# The following must have non-empty value if OpenMP is required
OMP :=

# The following must have non-empty value if Debugging compiler options are required
DEBUG :=

# Dependency Generator
DEPGEN := fortdepend
DEPGEN_INSTALL_DOCS := [https://github.com/ZedThree/fort_depend.py]

# Compiler Flags
FF += -O3 -march=native
ifeq ($(FC), gfortran)
  ifdef OMP
    FF += -fopenmp
  endif
else ifeq ($(FC), ifort)
  # FF += -heap-arrays
  ifdef DEBUG
  FF += -fp-stack-check -g -traceback -check bounds
  endif
  ifdef OMP
    FF += -qopenmp
  endif
endif

# Linker and linker Flags
LD := $(FC)
ifeq ($(LD), gfortran)
  LF += -fopenmp
else ifeq ($(LD), ifort)
  LF += -qopenmp
endif

# List of all source files
SRC_DIR := src
SRCS := $(notdir $(wildcard $(SRC_DIR)/*.f90))

# List of all object files
BUILD_DIR := build
OBJS := $(addprefix $(BUILD_DIR)/, $(addsuffix .o, $(basename $(SRCS))))

# Include path (to be searched by compiler for *.mod files)
IP := $(BUILD_DIR)
ifeq ($(FC), gfortran)
  FF += -J $(IP)
else ifeq ($(FC), ifort)
  FF += -module $(IP)
endif 

# Target executable(s)
EXECS := ccd_run ccd_cpt_to_xy ccd_rinit ccd_traj_to_legacy
EXECS := $(addprefix $(BUILD_DIR)/, $(EXECS))

# List of all executable scripts
SCRIPT_DIR := scripts
SCRIPTS := $(wildcard $(SCRIPT_DIR)/*)

# Dependency file to be generated using `fortdepend`
DEPFILE := .dependencies

# Intrinsic modules in standard Fortran for `fortdepend` to ignore
IMODS := omp_lib omp_lib_kinds iso_fortran_env ieee_arithmetic ieee_exceptions ieee_features iso_c_binding

# Font colors to be used by `echo`
RED='\e[1;31m'
GREEN='\e[1;32m'
BLUE='\e[1;34m'
NOCOLOR='\e[0m'

# Where to seek prerequisites
VPATH := $(SRC_DIR)

# System path where executables would be installed
INSTALL_PATH := /usr/local/bin

# Shell which runs the recipes
SHELL := bash

.PHONY: all clean rebuild install uninstall $(DEPGEN)

all: $(EXECS)
	@echo -e \\n$(GREEN)"make: Success"$(NOCOLOR)

$(EXECS): % : %.o $(filter-out $(BUILD_DIR)/ccd_%.o, $(OBJS))
	$(LD) $(LF) -o $@ $^
	@echo -e $(BLUE)"make: Built $@"$(NOCOLOR)

$(OBJS): $(BUILD_DIR)/%.o : %.f90
	$(FC) -c $(FF) -o $@ $<

# Rebuild all object files when this Makefile or dependency changes
# Create build directory only if non-existent (implemented as "order-only prerequisite")
$(OBJS): $(MAKEFILE_LIST) $(DEPFILE) | $(BUILD_DIR)

$(BUILD_DIR):
	mkdir $@

# Generate fresh dependency file whenever the codebase (sources) is modified or this Makefile changes
$(DEPFILE): $(SRCS) $(MAKEFILE_LIST) | $(DEPGEN)
	@echo 'Generating dependencies:'
	$(DEPGEN) --files $(addprefix $(SRC_DIR)/,$(SRCS)) --build $(BUILD_DIR) --ignore-modules $(IMODS) --output $(DEPFILE) --overwrite

# Define dependencies between object files
# Note: In fortran, object files are interdependent through .mod files
# Note: .mod file is generated only when module code is compiled into object file

include $(DEPFILE)

$(DEPGEN):
	@which $@ > /dev/null || { echo -e $(RED)"Make: Please install $@ first. $(DEPGEN_INSTALL_DOCS)"$(NOCOLOR) && false;}

clean: 
	rm -rf $(BUILD_DIR)
	rm -f $(DEPFILE)

rebuild: clean all

install:
	sudo install -t $(INSTALL_PATH) $(EXECS) $(SCRIPTS)

uninstall:
	sudo rm $(addprefix $(INSTALL_PATH)/, $(notdir $(EXECS)))
	sudo rm $(addprefix $(INSTALL_PATH)/, $(notdir $(SCRIPTS)))
