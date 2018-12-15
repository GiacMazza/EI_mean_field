#COMPILER (PARALLEL)
FC=mpif90
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

#EXE=q_dep_EI
EXE=mf_HKw90

HERE=`pwd`
DIR=${HERE}/drivers
DIREXE=${HOME}/.project_bin

.SUFFIXES: .f90

#
OBJS=VARS_GLOBAL.o HF.o 

LIBDIR=$(HOME)/opt_local
#LIBDIR=/opt/

#INCARGS=$(shell pkg-config --cflags --libs dmft_tools)
#INCARGS= -L/home/mazza/opt/dmft_tools/1.2.1/gnu/lib -ldmft_tools -I/home/mazza/opt/dmft_tools/1.2.1/gnu/include 
# INCARGS=-I$(LIBDIR)/update/SciFortran/gnu/include -L$(LIBDIR)/update/SciFortran/gnu/lib  -lscifor
# INCARGS+=$(shell pkg-config --cflags --libs dmft_tools)
#ARGS=-lscifor  

LIBARGS=$(shell pkg-config --libs   dmft_tools scifor)
INCARGS=$(shell pkg-config --cflags dmft_tools scifor)

#$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)

# Makefile.inc:INCARGS=$(shell pkg-config --cflags dmft_tools scifor)
# Makefile.inc:	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
# Makefile.inc:	$(FC) $(FPPFLAG) $(FFLAG) $(INCARGS) -c $<


FFLAG += -ffree-line-length-none  $(INCARGS) 


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
	$(FC)  $(FFLAG) $(INCARGS) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(LIBARGS)
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

