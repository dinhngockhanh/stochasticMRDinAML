.SUFFIXES: .f .f90 .o .mod
%.o : %.mod
FC		= gfortran
FFLAGS	= -mcmodel=medium -fdefault-integer-8 -O3 -fbacktrace

# FC		= ifort
#	-lacml
#	-fopenmp
# FFLAGS	= -mcmodel=large -O3 -heap-arrays -traceback

OBJLIBS 	= clock.o big_integer_module.o \
		  blas.o lapack.o slatec.o minpack.o \
		  expokit.o ode.o praxis.o nelmin.o randgen.o

main_OBJS  	= $(OBJLIBS) module_patient.o module_grand.o

program_fitting:			program_fitting.f $(main_OBJS)
	$(FC) $(FFLAGS) -o program_fitting program_fitting.f $(main_OBJS)

program_fitting_sensitive:	program_fitting_sensitive.f $(main_OBJS)
	$(FC) $(FFLAGS) -o program_fitting_sensitive program_fitting_sensitive.f $(main_OBJS)

program_stochastic:		program_stochastic.f $(main_OBJS)
	$(FC) $(FFLAGS) -o program_stochastic program_stochastic.f $(main_OBJS)

all: program_fitting program_fitting_sensitive program_stochastic

.f.o:;  $(FC) $(FFLAGS) -c $<
.f90.mod:;  $(FC) $(FFLAGS) -c $<
.f90.o:;  $(FC) $(FFLAGS) -c $<
