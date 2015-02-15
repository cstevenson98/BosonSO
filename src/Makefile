DEBUG=-O3

WARN=-Wall -Wno-unused -Wno-tabs

NAGDIR=-L/home/cstev/Libraries/fll6a24dfl/lib
NAGINC=-I/home/cstev/Libraries/fll6a24dfl/nag_interface_blocks
NAGLIB=-lnag_nag

LDFLAGS= $(DEBUG) $(NAGDIR)
FFLAGS= $(NAGINC)  $(DEBUG) $(WARN) 
LOADLIBES=$(NAGLIB)

FCOM=gfortran
LINK=gfortran

PROGRAMS=main 

all: $(PROGRAMS)
	mv *.o *.mod ../obj
	mv main ../
	@echo "DONE"

OBJS_MAIN=utility.o bosonsetup.o eom_so.o evolution.o main.o
main: $(OBJS_MAIN)
	$(LINK) $(LDFLAGS) $(OBJS_MAIN) -o $@ $(LOADLIBES) 	

ft=ft-data.o
ft: $(ft)
	$(LD) $(LDFLAGS) $(ft) -o $@ $(LOADLIBES) 	

%.o: %.f95
	$(FCOM) $(FFLAGS) -c $*.f95

clean:
	rm -f $(OBJS_MAIN) $(ODIR)*.mod $(PROGRAMS)

realclean: clean
	rm -f *~
