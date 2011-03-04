EXE = ./a.out

#MACRO = -DZFITTER# use the effective couplings, etc. computed by ZFITTER

SYS := $(shell uname -s)
ifeq ($(SYS),Darwin)
CC = /opt/local/bin/g++
PLATFORM = GNU-MacOSX
CFLAGS = -lgfortran
else
CC = g++
PLATFORM = GNU-Linux-x86
CFLAGS = -lgfortran
endif

ifeq ($(strip $(ROOTSYS)),)
$(error "Please set ROOTSYS.")
else
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
#ROOTLIBS += -lMinuit
#ROOTGLIBS = $(shell root-config --glibs) # for GUI
endif

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
EWPHYSICS_INC = -I../src
EWPHYSICS_LIB = -L../$(LIB_DIR)

INC_DIRS = $(ROOTCFLAGS) $(GSLPP_INC) $(UTILS_INC) $(SM_INC) $(ZFITTER6_43_INC) \
$(ZFITTER_INC) $(EWPHYSICS_INC)
LIB_DIRS = $(GSLPP_LIB) $(UTILS_LIB) $(SM_LIB) $(ZFITTER_LIB) $(EWPHYSICS_LIB)
LIBS = -lewphysics -lzfitter
ifeq ($(SYS),Darwin)
LIB_DIRS += $(ZFITTER6_43_LIB)
LIBS += -lzfitter6_43
else
#-- in case the commands ar and ranlib don't produce a correct library --
LIBS += ../../ZFitter6_43/src/*.o 
endif
LIBS += -lstandardmodel -lutils -lgslpp -lgslcblas -lgsl $(ROOTLIBS)

all: $(EXE)

$(EXE): EWphysicsTest.cpp Makefile
	$(CC) $(MACRO) -o $@ $(INC_DIRS) $< $(CFLAGS) $(LIB_DIRS) $(LIBS)

clean:
	rm -f $(EXE)