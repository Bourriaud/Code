OBJECTS = main.o inout.o efficiency.o

run: $(OBJECTS)
	gfortran -o run $(OBJECTS)

main.o : main.f90 inout.mod efficiency.mod
	gfortran -c main.f90

inout.mod : inout.o

inout.o : inout.f90
	gfortran -c inout.f90

efficiency.mod : efficiency.o

efficiency.o : efficiency.f90
	gfortran -c efficiency.f90

clean :
	rm -f *.o *.mod *~
