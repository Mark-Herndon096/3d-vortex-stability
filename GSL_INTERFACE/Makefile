.SUFFIXES:

# Object list
OBJECTS = special_function_wrapper.o \
	  special_function_interface.o  \

CFLAGS = -O2

MODDIR = ../.mod 
FCFLAGS = -module $(MODDIR) 

%.o: %.c; $(CC) $(CFLAGS) -c -o $@ $< $(LIBRARIES)
%.o: %.f90; $(FC) $(COMPFLAGS) $(FCFLAGS) -c -o $@ $< $(LIBRARIES)

gsl_objs: $(OBJECTS)

clean:
	rm -rf *.o *.mod $(EXEC_NAME)

# Object Dependencies
special_function_wrapper.o: special_function_wrapper.c
special_function_interface.o: special_function_interface.f90 special_function_wrapper.o

