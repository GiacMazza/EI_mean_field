#TO BE CHANGED BY USER:
#$ DRIVER NAME without .f90 extension
#$ COMPILER: supported compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default if $HOME/.bin in the PATH)

#EXE=q_dep_EI
#EXE=mf_HKw90
#EXE=3b_cfr_kaneko
#EXE=3b_simple
#EXE=3b_full
#EXE=toy_tns
#EXE=toy_tns_ei
#EXE=tns_symmetry_breaking
EXE=two_chains_tns
EXE=two_chains_dynamics
#EXE=two_chains_energy_fixed_phase
#EXE=tns_6bands
#EXE=two_chains_restric_hf
#EXE=fulltwo_chains_tns
#EXE=two_chains_cmplx
#EXE=sweep_V_two_chains_tns


FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin
DIR=drivers

define colorecho	
	@tput setaf 6
	@tput bold
	@echo $1
	@tput sgr0
endef


#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC:=$(shell pkg-config --cflags dmft_tools dmft_ed  scifor)
GLOB_LIB:=$(shell pkg-config --libs dmft_ed dmft_tools scifor)



ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif
ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising -Wuninitialized -fbounds-check  -Waliasing -Wall -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI	
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc


OBJS=VARS_GLOBAL.o R_HF.o


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

.f90.o:
	$(FC) $(FLAG) $(GLOB_INC) -c $<


all: FLAG:=${FFLAG} ${FPPMPI}
all:	$(OBJS)
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

debug: FLAG:=${DFLAG} ${FPPMPI}
debug: $(OBJS)
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"



serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

serial_debug: FLAG:=${DFLAG} ${FPPSERIAL}
serial_debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXE)



#########################################################################
