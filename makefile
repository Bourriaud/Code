OBJECTS = main.o inout.o efficiency.o types.o phys.o

run: $(OBJECTS)
	gfortran -o run $(OBJECTS)

main.o : main.f90 inout.mod efficiency.mod types.mod phys.mod
	gfortran -c main.f90

inout.mod : inout.o

inout.o : inout.f90 types.mod
	gfortran -c inout.f90

efficiency.mod : efficiency.o

efficiency.o : efficiency.f90 types.mod
	gfortran -c efficiency.f90

types.mod : types.o

types.o : types.f90
	gfortran -c types.f90

phys.mod : phys.o

phys.o : phys.f90
	gfortran -c phys.f90

clean :
	rm -f *.o *.mod *~
