# Collective Cell Dynamics (ccd)

FC=ifort
FF= -O3
EXE=ccd.exe

.PHONY : all clean

all : $(EXE) clean

$(EXE) : main.f90
	@$(FC) $(FF) -o $(EXE) main.f90

clean :
	@rm -f *.o *.mod
