PROG = run
SRCS = main.f90 constant.f90 inout.f90 efficiency.f90 types.f90 phys.f90 FV.f90 time.f90 reconstruction.f90
OBJS = main.o constant.o inout.o efficiency.o types.o phys.o FV.o time.o reconstruction.o

F90 = gfortran
F90FLAGS=-O0 -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic
F90FLAGS2=-O0 -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic -g -pg
F90FLAGS3=-O2 -march=native

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

main.o: constant.o inout.o efficiency.o types.o phys.o FV.o time.o reconstruction.o

inout.o: constant.o types.o

efficiency.o: constant.o types.o

types.o: constant.o

phys.o: constant.o

FV.o: constant.o types.o phys.o inout.o efficiency.o

time.o: constant.o types.o FV.o reconstruction.o

reconstruction.o: constant.o types.o

test: test.o constant.o reconstruction.o types.o inout.o
	$(F90) -o test test.o reconstruction.o types.o inout.o

test.o : constant.o reconstruction.o types.o inout.o
