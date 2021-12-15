# Makefile for primary src
.SUFFIXES:

# Object list
OBJECTS = mod_global.o               \
	  mod_numerical_routines.o   \
	  mod_io.o                   \
	  main.o

%.o: %.f90; $(FC) $(COMPFLAGS) -c -o $@ $< $(LIBRARIES)

src_objs: $(OBJECTS)

clean:
	rm -rf *.o *.mod $(EXEC_NAME)

## Object dependencies
#mod_global.o: mod_global.f90
#mod_numerical_routines.o: mod_numerical_routines.f90
#mod_file_io.o: mod_file_io.f90 mod_global.o
#main.o: main.f90 mod_file_io.o mod_global.o mod_numerical_routines.o
# Makefile:1 ends here
