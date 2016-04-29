# Program for 3d electrodeposition simulation-Salida NETCDF
# Esteban Mocskos - 2002 -2007

C = gcc
#Version gcc
CPP = g++
#CPP = icc

#Para Itanium(gcc):
#CFLAGS_ORI_prev= -O3 -Wall -Winline -D$(PROGRAM) -I /include -I $(STRFUNC)

#Para Itanium(icc):
#CFLAGS= -mcpu=itanium2 -fast -gcc-name=gcc-4.1 
#CFLAGS= -IPF_fp_speculation fast  -IPF_fma -mcpu=itanium2 -ftz -tpp2 -D$(PROGRAM) -I /include -I $(NETCDFINC)  -I $(STRFUNC)

#Para profilling
#CFLAGS_ORI_prev= -Wall -Winline -D$(PROGRAM) -pg

#Pentium 4(gcc)
# CFLAGS_ORI_prev= -march=pentium4 -O3 -mfpmath=sse -msse2 -mmmx -Wall -fomit-frame-pointer -pipe -DMPICH_IGNORE_CXX_SEEK
# -DDIFFERENTIATED_DOMAINS
# -DSIN_CONVECCION

#Pentium 4(icc)
#CFLAGS= -tpp7 -xW -axW -march=pentium4 -D$(PROGRAM) -I/usr/local/include .0
#CFLAGS= -tpp7 -mcpu=pentium4 -D$(PROGRAM) -I/usr/local/include -ansi

#Para las opterons
#CFLAGS= -march=opteron -O3 -mfpmath=sse -msse2 -mmmx -Wall -fomit-frame-pointer -pipe -DMPICH_IGNORE_CXX_SEEK

#Para las Xeon64EMT(gcc)
# CFLAGS_ORI_prev= -O3 -march=nocona -mfpmath=sse -msse2 -mmmx -Wall -fomit-frame-pointer -pipe -DMPICH_IGNORE_CXX_SEEK
#CFLAGS_ORI_prev= -fast -parallel -msse2 -Wall 
# -DDIFFERENTIATED_DOMAINS
# -DSIN_CONVECCION

# Si se quiere hacer el profiling
# CFLAGS_ORI= -O0 -march=nocona -mfpmath=sse -msse2 -mmmx -Wall -pipe -pg 


# Para tener soporte de X, dejar esta variable declarada
# WITH_WINDOWS_SUPP = "dummy"

#Normal.
CFLAGS= -O3 -Wall
#Sin optimizacion.
# CFLAGS= -Wall -g 

#ifdef WITH_WINDOWS_SUPP
#CFLAGS_ORI = $(CFLAGS_ORI_prev) -DWITH_WINDOWS_SUPP
#else
#CFLAGS_ORI = $(CFLAGS_ORI_prev)
#endif

#ifdef PROGRAM
#CFLAGS= $(CFLAGS_ORI) -D$(PROGRAM)
#endif

# Normal
CLINKFLAGS = -gcc-name=gcc-4.4 -gxx-name=g++-4.4 

# Para Profiling
# CLINKFLAGS = -pg

LIBS = -lgsl -lgslcblas -lstdc++

NOMBRE = main

# -lstrfunc 
#OBJS= main.o Entorno.o Potential.o ElectricField.o Specie.o 
OBJS= main.o Entorno.o ElectricField.o Potential.o Specie.o 
#OBJS= main.o Entorno.o Potential.o Specie.o Optimizacion.o
.SUFFIXES: .cpp .o

.cpp.o:
#	echo "Compiling $@";$(CPP) $(CFLAGS) -c $<
	$(CPP) $(CFLAGS) -c $<

.c.o:
#	@echo "Compiling $@";$(C) $(CFLAGS) -c $<
	$(C) $(CFLAGS) -c $<

program: $(OBJS); 
	echo "Linking";${CPP} ${CLINKFLAGS} ${OBJS} ${LIBS} -o ${NOMBRE} 

clean:
	@echo "Removing Object files.";rm -f *.o; 

all:
	@(make clean; make;)

default:
	@(echo "making $(PROGRAM)"; export PROGRAM=054d1.4D1.2.324.80; $(MAKE) program;)
