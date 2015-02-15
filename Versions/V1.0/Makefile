DEBUG=-O3

WARN=-Wall -Wno-unused -Wno-tabs

NAGDIR=-L/home/conchops/Libraries/fll6a24dfl/lib
NAGINC=-I/home/conchops/Libraries/fll6a24dfl/nag_interface_blocks
NAGLIB=-lnag_nag

LDFLAGS= $(DEBUG) $(NAGDIR)
FFLAGS= $(NAGINC)  $(DEBUG) $(WARN) 
LOADLIBES=$(NAGLIB)

FC=gfortran
LD=gfortran

PROGRAMS=main 

all: $(PROGRAMS)
	@echo "DONE"
	
main=utility.o boseparams.o state_control.o main.o
main: $(main)
	$(LD) $(LDFLAGS) $(main) -o $@ $(LOADLIBES) 	
	
TEMPLATE=
TEMPLATE: $(TEMPLATE)
	$(LD) $(LDFLAGS) $(TEMPLATE) -o $@ $(LOADLIBES) 	

%.o: %.f95
	$(FC) $(FFLAGS) -c $*.f95

clean:
	rm -f *.o *.mod $(PROGRAMS)

realclean: clean
	rm -f *~
