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
FC=gfortran-mp-4.5
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/ZFitter.o \
	${OBJECTDIR}/src/ZFitterObservables.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1 \
	${TESTDIR}/TestFiles/f2

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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a

${OBJECTDIR}/src/ZFitter.o: src/ZFitter.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ZFitter.o src/ZFitter.cpp

${OBJECTDIR}/src/ZFitterObservables.o: src/ZFitterObservables.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ZFitterObservables.o src/ZFitterObservables.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/ZFitterTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}  -lgsl -lgslcblas -L/usr/lib/root -L/usr/local/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} -lgfortran ../ZFitter6_43/dist/Debug/GNU-MacOSX/libzfitter6_43.a ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a ../StandardModel/dist/Debug/GNU-MacOSX/libstandardmodel.a 

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/ZFitterTest-2.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}  -lgsl -lgslcblas -L/usr/lib/root -L/usr/local/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} -lgfortran ../ZFitter6_43/dist/Debug/GNU-MacOSX/libzfitter6_43.a ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a ../StandardModel/dist/Debug/GNU-MacOSX/libstandardmodel.a 


${TESTDIR}/tests/ZFitterTest.o: tests/ZFitterTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -Isrc -I../gslpp/src -I../StandardModel/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/ZFitterTest.o tests/ZFitterTest.cpp


${TESTDIR}/tests/ZFitterTest-2.o: tests/ZFitterTest-2.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -Isrc -I../gslpp/src -I../StandardModel/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/ZFitterTest-2.o tests/ZFitterTest-2.cpp


${OBJECTDIR}/src/ZFitter_nomain.o: ${OBJECTDIR}/src/ZFitter.o src/ZFitter.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ZFitter.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ZFitter_nomain.o src/ZFitter.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ZFitter.o ${OBJECTDIR}/src/ZFitter_nomain.o;\
	fi

${OBJECTDIR}/src/ZFitterObservables_nomain.o: ${OBJECTDIR}/src/ZFitterObservables.o src/ZFitterObservables.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ZFitterObservables.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ZFitterObservables_nomain.o src/ZFitterObservables.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ZFitterObservables.o ${OBJECTDIR}/src/ZFitterObservables_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libzfitter.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
