#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/EWphysics.o

# Test Directory
TESTDIR=build/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-Release.mk dist/Release/GNU-MacOSX/libewphysics.a

dist/Release/GNU-MacOSX/libewphysics.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-MacOSX
	${RM} dist/Release/GNU-MacOSX/libewphysics.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewphysics.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/GNU-MacOSX/libewphysics.a

${OBJECTDIR}/src/EWphysics.o: src/EWphysics.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWphysics.o src/EWphysics.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}

${TESTDIR}/tests/EWphysicsTest.o: tests/EWphysicsTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -MMD -MP -MF $@.d -o ${TESTDIR}/tests/EWphysicsTest.o tests/EWphysicsTest.cpp


${OBJECTDIR}/src/EWphysics_nomain.o: ${OBJECTDIR}/src/EWphysics.o src/EWphysics.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWphysics.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWphysics_nomain.o src/EWphysics.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWphysics.o ${OBJECTDIR}/src/EWphysics_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/Release
	${RM} dist/Release/GNU-MacOSX/libewphysics.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
