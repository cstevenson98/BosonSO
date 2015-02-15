HOSTNAME := $(shell hostname)

include Makefile.$(HOSTNAME)

DEBUG=-fbacktrace -ggdb -fbounds-check 
#DEBUG=-O3

WARN=-Wall -Wno-unused

LDFLAGS= $(DEBUG) $(NAGDIR)  $(LAPACKDIR)
FFLAGS= $(NAGINC)  $(DEBUG) $(WARN) 
LOADLIBES=$(NAGLIB)


FC=gfortran
LD=gfortran

PROGRAMS=gh-dynamics rs-dynamics

all: $(PROGRAMS)
	@echo "DONE"

TEMPLATE=
TEMPLATE: $(TEMPLATE)
	$(LD) $(LDFLAGS) $(TEMPLATE) -o $@ $(LOADLIBES) 

COMMON=utility.o run-diffeqn.o params.o hermite.o

gh=$(COMMON) gh-overlaps.o gh-eom.o gh-dynamics.o 
gh-dynamics: $(gh)
	$(LD) $(LDFLAGS) $(gh) -o $@ $(LOADLIBES) 

rs=$(COMMON) rs-basis.o rs-eom.o rs-dynamics.o 
rs-dynamics: $(rs)
	$(LD) $(LDFLAGS) $(rs) -o $@ $(LOADLIBES) 




%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

%.mod: %.o


%.f90:
	ln -s ~/common-code/$@ ./

clean:
	rm -f *.o *.mod $(PROGRAMS)

realclean: clean
	rm -f *~
