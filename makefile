OBJECTS = main.o constant.o inout.o efficiency.o types.o phys.o FV.o time.o reconstruction.o

run: $(OBJECTS)
	gfortran -o run $(OBJECTS)

main.o : main.f90 constant.mod inout.mod efficiency.mod types.mod phys.mod FV.mod time.mod reconstruction.mod
	gfortran -c main.f90

constant.mod : constant.o

constant.o : constant.f90
	gfortran -c constant.f90

inout.mod : inout.o

inout.o : inout.f90 constant.mod types.mod
	gfortran -c inout.f90

efficiency.mod : efficiency.o

efficiency.o : efficiency.f90 constant.mod types.mod
	gfortran -c efficiency.f90

types.mod : types.o

types.o : types.f90 constant.mod
	gfortran -c types.f90

phys.mod : phys.o

phys.o : phys.f90 constant.mod
	gfortran -c phys.f90

FV.mod : FV.o

FV.o : FV.f90 constant.mod types.mod phys.mod inout.mod efficiency.mod
	gfortran -c FV.f90

time.mod : time.o

time.o : time.f90 constant.mod types.mod FV.mod reconstruction.mod
	gfortran -c time.f90

reconstruction.mod : reconstruction.o

reconstruction.o : reconstruction.f90 constant.mod types.mod
	gfortran -c reconstruction.f90

clean :
	rm -f *.o *.mod *~

test: test.o constant.o reconstruction.o types.o inout.o
	gfortran -o test test.o reconstruction.o types.o inout.o

test.o : test.f90 constant.mod reconstruction.mod types.mod inout.mod
	gfortran -c test.f90
