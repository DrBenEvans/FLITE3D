#Compilers and flags
#for Linux FC4

MF	= ifort 
RM	= /bin/rm

#FFLAGS	= -r8 -g -C  -traceback
FFLAGS	= -r8 -O2   -traceback
FFLAGS90 = 

#Fortran Object Files

FOBJS   = \
          gen3d_tools.o  \
          gen3d_struct.o \
          gen3d_elements.o \
          gen3d_boundary.o \
          gen3d_grid.o \
          gen3d_generator.o \
          gen3d_main.o \

OBJS	= $(FOBJS)


#Compile and link

default:
	make all

all:  $(OBJS)
	$(MF) $(FFLAGS) -o Gen3d  $(OBJS) -L/eng/cvcluster/egevansbj/codes/prepro/Metis -lmetis
		
clean:
	$(RM) -f $(OBJS) Gen3d *.mod

gen3d_tools.o:gen3d_tools.f90
	$(MF) $(FFLAGS) -c gen3d_tools.f90
gen3d_struct.o:gen3d_struct.f90
	$(MF) $(FFLAGS) -c gen3d_struct.f90
gen3d_elements.o:gen3d_elements.f90
	$(MF) $(FFLAGS) -c gen3d_elements.f90
gen3d_boundary.o:gen3d_boundary.f90
	$(MF) $(FFLAGS) -c gen3d_boundary.f90
gen3d_grid.o:gen3d_grid.f90
	$(MF) $(FFLAGS) -c gen3d_grid.f90
gen3d_generator.o:gen3d_generator.f90
	$(MF) $(FFLAGS) -c gen3d_generator.f90
gen3d_main.o:gen3d_main.f90
	$(MF) $(FFLAGS) -c gen3d_main.f90

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

