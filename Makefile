#COMPILER (PARALLEL)
FC=mpif90
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

EXE=q_dep_EI

HERE=`pwd`
DIR=${HERE}/drivers
DIREXE=${HOME}/.project_bin

.SUFFIXES: .f90

#
OBJS=VARS_GLOBAL.o 

#GALLIBDIR  = /home/mazza/opt_local/galahad/objects/pc64.lnx.gfo/double
# GALLIBDIR  = /opt/galahad/objects/pc64.lnx.gfo/double
# GALLIBS1   = -lgalahad -lgalahad_hsl 
# GALLIBS2   = -lgalahad_metis 

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#FFLAG +=-fpp -D_$(FPP)
LIBDIR=$(HOME)/opt_local
#LIBDIR=/opt/

INCARGS=-I$(LIBDIR)/SciFortran/gnu/include -L$(LIBDIR)/SciFortran/gnu/lib 
INCARGS+=-I$(LIBDIR)/DMFTtools/gnu/include -L$(LIBDIR)/DMFTtools/gnu/lib 
FFLAG += -ffree-line-length-none -cpp $(INCARGS)
#FFLAG+=-O0 -p -g -Wall -fbounds-check -fbacktrace -Wuninitialized

#ARGS=  $(GALLIBS1) $(GALLIBS2) -I/home/mazza/opt_local/galahad/modules/pc64.lnx.gfo/double -lscifor $(MKLARGS) -lminpack -larpack -lparpack   
ARGS=-ldmftt -lscifor  -lfftpack -lminpack  -llapack -lblas -larpack 


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
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " !+--------------------------------------+! "
	@echo " .................. done .................. "
	@echo " !+--------------------------------------+! "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

ed_solver:
	@make -C ED_SOLVER/

.f90.o:	
	$(FC) $(FFLAG)  -c $< 

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

