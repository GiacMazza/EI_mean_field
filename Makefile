#COMPILER (PARALLEL)
FC=mpif90
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

#EXE=q_dep_EI
#EXE=mf_HKw90
#EXE=3b_cfr_kaneko
#EXE=3b_simple
#EXE=3b_full
#EXE=toy_tns
#EXE=toy_tns_ei
#EXE=tns_symmetry_breaking
#EXE=two_chains_tns
EXE=tns_6bands
#EXE=two_chains_restric_hf
#EXE=fulltwo_chains_tns
#EXE=two_chains_cmplx
#EXE=sweep_V_two_chains_tns


DIR=./drivers
DIREXE=${HOME}/.bin

.SUFFIXES: .f90

#
#OBJS=VARS_GLOBAL.o HF.o  R_HF.o
OBJS=VARS_GLOBAL.o R_HF.o

LIBDIR=$(HOME)/gm_opt
#LIBDIR=/opt/



INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
LIBARGS=$(shell pkg-config --libs   dmft_tools)
LIBARGS+=$(shell pkg-config --libs   scifor)


scifor_DIR=$(HOME)/gm_opt/scifor/gnu/4.7.4
dmft_tools_DIR=$(HOME)/gm_opt/dmft_tools/gnu/2.3.0-10-g6ae56d5
INCARGS =  -I$(dmft_tools_DIR)/include -I$(scifor_DIR)/include
LIBARGS =  -L$(dmft_tools_DIR)/lib   -L$(scifor_DIR)/lib -L/usr/lib64  -ldmft_tools  -lscifor -llapack -lblas 


#$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)

# Makefile.inc:INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
# Makefile.inc:	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
# Makefile.inc:	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) -c $<


FFLAG += -ffree-line-length-none  $(INCARGS) 
FFLAG+=-O0 -p -g -Wall -fbounds-check -fbacktrace -Wuninitialized

BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
ifeq ($(BRANCH),_master)
BRANCH=
endif

all:compile


lib: ed_solver

compile: version $(OBJS)
	@echo " !+------------------------------------------------- "
	@echo " ..................... compile ..................... "
	@echo " !+------------------------------------------------- "
	$(FC)  $(FFLAG)  $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
	@echo " !+--------------------------------------+! "
	@echo " .................. done .................. "
	@echo " !+--------------------------------------+! "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

ed_solver:
	@make -C ED_SOLVER/

.f90.o:	
	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) -c $<	

completion:
	sf_lib_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

all_clean: clean
	@make -C ED_SOLVER/ clean

version:
	@echo $(VER)

