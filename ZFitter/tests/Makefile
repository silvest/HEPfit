EXE = ./a.out

SYS := $(shell uname -s)
ifeq ($(SYS),Darwin)
PLATFORM = GNU-MacOSX
CFLAGS = -lgfortran
else
PLATFORM = GNU-Linux-x86
CFLAGS = -lgfortran
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
ZFITTER_INC = -I../src
ZFITTER_LIB = -L../$(LIB_DIR) 

INC_DIRS = $(GSLPP_INC) $(UTILS_INC) $(SM_INC) $(ZFITTER6_43_INC) $(ZFITTER_INC)
LIB_DIRS = $(GSLPP_LIB) $(UTILS_LIB) $(SM_LIB) $(ZFITTER_LIB)
ifeq ($(SYS),Darwin)
LIB_DIRS += $(ZFITTER6_43_LIB)
LIBS = -lzfitter6_43
else
LIBS = ../../ZFitter6_43/src/*.o # the commands ar and ranlib don't produce a correct library.
endif
LIBS += -lzfitter -lstandardmodel -lutils -lgslpp -lgslcblas -lgsl

all: $(EXE)

$(EXE): ZFitterTest.cpp Makefile
	g++ -o $@ $(INC_DIRS) $< $(CFLAGS) $(LIB_DIRS) $(LIBS)

clean:
	rm -f $(EXE)