#Compilers and flags
#for Linux FC4

MF	= mpiifort 
RM	= /bin/rm

FFLAGS	= -r8 -cpp -O2 
#FFLAGS	= -r8 -cpp -g -traceback  -I$(MPI_HOME)/include
FFLAGS90 = 
FFLAGSMP = -DPARALLEL -DMASTER 

#Fortran Object Files

FOBJS   = \
          mg3d_mpilib.o \
          mgnsg3d_def_par.o \
          mgnsg3d_mg_par.o \
          mgnsg3d_main_par.o

OBJS	= $(FOBJS)


#Compile and link

default:
	make all

all:  $(OBJS)
	$(MF) $(FFLAGS) -o UnsMgnsg3d $(OBJS) 
		
clean:
	$(RM) -f $(OBJS) UnsMgnsg3d *.mod

mg3d_mpilib.o:mg3d_mpilib.F90
	$(MF) $(FFLAGS) -c mg3d_mpilib.F90
mgnsg3d_def_par.o:mgnsg3d_def_par.F90
	$(MF) $(FFLAGS) $(FFLAGS90) $(FFLAGSMP) -c mgnsg3d_def_par.F90
mgnsg3d_mg_par.o:mgnsg3d_mg_par.F90
	$(MF) $(FFLAGS) $(FFLAGS90) $(FFLAGSMP) -c mgnsg3d_mg_par.F90
mgnsg3d_main_par.o:mgnsg3d_main_par.F90
	$(MF) $(FFLAGS) $(FFLAGSMP)-c mgnsg3d_main_par.F90

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

