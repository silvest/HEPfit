EXE = ./a.out

#MACRO = -DZFITTER# use the effective couplings, etc. computed by ZFITTER

SYS := $(shell uname -s)
ifeq ($(SYS),Darwin)
CC = /usr/bin/g++
#CC = /opt/local/bin/g++
#CC = /Developer/usr/bin/g++
PLATFORM = GNU-MacOSX
CFLAGS = -lgfortran
#BOOSTLIBRARY_INC = -L/opt/local/include
#BOOSTLIBRARY_LIB = -L/opt/local/lib
else
CC = g++
PLATFORM = GNU-Linux-x86
CFLAGS =  -g -Wall -lgfortran
#BOOSTLIBRARY_INC = 
#BOOSTLIBRARY_LIB = 
endif

#ifeq ($(strip $(ROOTSYS)),)
#$(error "Please set ROOTSYS.")
#else
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
#ROOTLIBS += -lMinuit
#ROOTGLIBS = $(shell root-config --glibs) # for GUI
#endif

LIB_DIR = dist/Debug/$(PLATFORM)
#LIB_DIR = dist/Release/$(PLATFORM)

GSLPP_INC = -I../../gslpp/src
GSLPP_LIB = -L../../gslpp/$(LIB_DIR) 
UTILS_INC = -I../../Utils/src
UTILS_LIB = -L../../Utils/$(LIB_DIR) 
SM_INC = -I../../StandardModel/src
SM_LIB = -L../../StandardModel/$(LIB_DIR) 
ZFITTER6_43_INC = -I../../ZFitter6_43/src
ZFITTER6_43_LIB = -L../../ZFitter6_43/lib 
ZFITTER_INC = -I../../ZFitter/src
ZFITTER_LIB = -L../../ZFitter/$(LIB_DIR)
LOOPFUNCTIONS_INC = -I../src
LOOPFUNCTIONS_LIB = -L../$(LIB_DIR)


INC_DIRS = $(GSLPP_INC) $(LOOPFUNCTIONS_INC) #$(BOOSTLIBRARY_INC)
LIB_DIRS = $(GSLPP_LIB) $(LOOPFUNCTIONS_LIB) #$(BOOSTLIBRARY_LIB)

LIBS = -lloopfunctions -lgslpp -lgslcblas -lgsl
#LIBS += -lboost_program_options

all: $(EXE)

$(EXE): LoopFunctionsTest.cpp Makefile
	$(CC) $(MACRO) -o $@ $(INC_DIRS) $< $(CFLAGS) $(LIB_DIRS) $(LIBS)

clean:
	rm -f $(EXE)
