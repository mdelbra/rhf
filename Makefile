# RHF Makefile 3 Jun 2014 
# by mdelbra

# Standard CC
CC=cc
CXX=c++

#By default no OMP
#OMP=0

# C optimization flags
COPT	= -O3 -Wall

# Linking flags
LDFLAGS   =  -lIlmImf -lHalf -lpng -pthread

ifeq (${OMP},1)
    COPT += -fopenmp
    LDFLAGS += -fopenmp
endif

# Lib PNG
PNG_INCLUDE=/usr/local/include
PNG_LIB=/usr/local/lib

# OpenEXR
EXR_INCLUDE1=/usr/local/include/OpenEXR
EXR_INCLUDE2=/usr/include/OpenEXR/
EXR_LIB=/usr/local/lib

# C++ optimization flags
CXXOPT	= $(COPT)

CFLAGS    = $(COPT) -I${EXR_INCLUDE1} -I${EXR_INCLUDE2} -I${PNG_INCLUDE}
CXXFLAGS  = $(CXXOPT) -I${EXR_INCLUDE1} -I${EXR_INCLUDE2} -I${PNG_INCLUDE}

EXEC      = rhf exrdiff exrcrop exrtopng

default: all

all: rhf exrdiff exrcrop exrtopng


# ------- C++ files -------
./io_png.o: ./io_png.c ./io_png.h
	$(CC) $(CFLAGS) -c ./io_png.c -o ./io_png.o

./exrcrop.o: ./exrcrop.cpp
	$(CXX) $(CXXFLAGS) -c ./exrcrop.cpp -o ./exrcrop.o

./exrtopng.o: ./exrtopng.cpp
	$(CXX) $(CXXFLAGS) -c ./exrtopng.cpp -o ./exrtopng.o

./exrdiff.o: ./exrdiff.cpp
	$(CXX) $(CXXFLAGS) -c ./exrdiff.cpp -o ./exrdiff.o

./io_exr.o: ./io_exr.cpp ./io_exr.h
	$(CXX) $(CXXFLAGS) -c ./io_exr.cpp -o ./io_exr.o

./libauxiliar.o: ./libauxiliar.cpp ./libauxiliar.h
	$(CXX) $(CXXFLAGS) -c ./libauxiliar.cpp -o ./libauxiliar.o

./libdenoising.o: ./libdenoising.cpp ./libdenoising.h
	$(CXX) $(CXXFLAGS) -c ./libdenoising.cpp -o ./libdenoising.o

./rhf.o: ./rhf.cpp
	$(CXX) $(CXXFLAGS) -c ./rhf.cpp -o ./rhf.o


# ------- Main -------
exrdiff: ./io_exr.o ./exrdiff.o
	$(CXX)  -L${EXR_LIB} -L${PNG_LIB}  ./io_exr.o  ./exrdiff.o  $(LDFLAGS) -o exrdiff

exrcrop:  ./io_exr.o ./exrcrop.o
	$(CXX)  -L${EXR_LIB} -L${PNG_LIB} ./io_exr.o  ./exrcrop.o  $(LDFLAGS) -o exrcrop


exrtopng:  ./io_exr.o ./io_png.o ./exrtopng.o 
	$(CXX)  -L${EXR_LIB} -L${PNG_LIB} ./io_exr.o  ./io_png.o  ./exrtopng.o  $(LDFLAGS)   -o exrtopng

rhf:  ./io_exr.o ./libauxiliar.o ./libdenoising.o ./rhf.o
	$(CXX)  -L${EXR_LIB} -L${PNG_LIB} ./io_exr.o  ./libauxiliar.o ./libdenoising.o ./rhf.o  $(LDFLAGS)   -o rhf



.PHONY: clean distclean

clean: 
	rm -f *.o

distclean: clean
	rm -f $(EXEC)

