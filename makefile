PROG = run
SRCS = main.f90 constant.f90 inout.f90 efficiency.f90 types.f90 phys.f90 FV.f90 time.f90 reconstruction.f90 limit.f90 ICBC.f90 AMR.f90 p4est_wrapper.c
OBJS = main.o constant.o inout.o efficiency.o types.o phys.o FV.o time.o reconstruction.o limit.o ICBC.o AMR.o p4est_wrapper.o

INC = -I/home/imb/abourriaud/Documents/Code/P4EST/local/include
LIBP = -L/home/imb/abourriaud/Documents/Code/P4EST/local/lib
LIBS = -lp4est -lsc -lmpi

CC = gcc
CCFLAGS = -O0

F90 = gfortran
F90FLAGS = -O0 -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic
F90FLAGS2 = -O0 -Wall -ffpe-trap=invalid,overflow -fcheck=all -pedantic -g -pg
F90FLAGS3 = -O2 -march=native

include /home/imb/abourriaud/Documents/Code/P4EST/Makefile.p4est.mk

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS) $(INC) $(LIBP) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod test.o test

.SUFFIXES: $(SUFFIXES) .c

.c.o:
	$(CC) $(CCFLAGS) -c $<  $(INC)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

main.o: constant.o inout.o efficiency.o types.o phys.o FV.o time.o reconstruction.o ICBC.o AMR.o

inout.o: constant.o types.o efficiency.o phys.o

efficiency.o: constant.o types.o phys.o

types.o: constant.o

phys.o: constant.o

FV.o: constant.o types.o phys.o inout.o efficiency.o reconstruction.o

time.o: constant.o types.o FV.o reconstruction.o limit.o

reconstruction.o: constant.o types.o efficiency.o

limit.o: constant.o types.o reconstruction.o

ICBC.o: constant.o types.o phys.o

AMR.o: constant.o types.o

test: test.o constant.o reconstruction.o types.o inout.o efficiency.o phys.o FV.o time.o limit.o ICBC.o p4est_wrapper.o
	$(F90) -o test constant.o test.o reconstruction.o types.o inout.o efficiency.o phys.o FV.o time.o limit.o ICBC.o AMR.o p4est_wrapper.o $(INC) $(LIBP) $(LIBS)

test.o : constant.o reconstruction.o types.o inout.o efficiency.o phys.o FV.o time.o limit.o ICBC.o AMR.o
