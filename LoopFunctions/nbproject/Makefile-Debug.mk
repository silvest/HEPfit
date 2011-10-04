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
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Polylogarithms.o \
	${OBJECTDIR}/src/PVfunctions.o \
	${OBJECTDIR}/src/ClausenFunctions.o \
	${OBJECTDIR}/src/BernoulliNumbers.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f3 \
	${TESTDIR}/TestFiles/f2

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-Wall
CXXFLAGS=-Wall

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a

${OBJECTDIR}/src/Polylogarithms.o: src/Polylogarithms.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Polylogarithms.o src/Polylogarithms.cpp

${OBJECTDIR}/src/PVfunctions.o: src/PVfunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/PVfunctions.o src/PVfunctions.cpp

${OBJECTDIR}/src/ClausenFunctions.o: src/ClausenFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ClausenFunctions.o src/ClausenFunctions.cpp

${OBJECTDIR}/src/BernoulliNumbers.o: src/BernoulliNumbers.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BernoulliNumbers.o src/BernoulliNumbers.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f3: ${TESTDIR}/tests/LFtestclass.o ${TESTDIR}/tests/LFtestrunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f3 $^ ${LDLIBSOPTIONS} -lcppunit -lgsl -lgslcblas ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a 

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/LoopFunctionsTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} -lgsl -lgslcblas ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a 


${TESTDIR}/tests/LFtestclass.o: tests/LFtestclass.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -I. -I. -Isrc -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/LFtestclass.o tests/LFtestclass.cpp


${TESTDIR}/tests/LFtestrunner.o: tests/LFtestrunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -I. -I. -Isrc -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/LFtestrunner.o tests/LFtestrunner.cpp


${TESTDIR}/tests/LoopFunctionsTest.o: tests/LoopFunctionsTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -I. -I. -I. -I. -I. -Isrc -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/LoopFunctionsTest.o tests/LoopFunctionsTest.cpp


${OBJECTDIR}/src/Polylogarithms_nomain.o: ${OBJECTDIR}/src/Polylogarithms.o src/Polylogarithms.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/Polylogarithms.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Polylogarithms_nomain.o src/Polylogarithms.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/Polylogarithms.o ${OBJECTDIR}/src/Polylogarithms_nomain.o;\
	fi

${OBJECTDIR}/src/PVfunctions_nomain.o: ${OBJECTDIR}/src/PVfunctions.o src/PVfunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/PVfunctions.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/PVfunctions_nomain.o src/PVfunctions.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/PVfunctions.o ${OBJECTDIR}/src/PVfunctions_nomain.o;\
	fi

${OBJECTDIR}/src/ClausenFunctions_nomain.o: ${OBJECTDIR}/src/ClausenFunctions.o src/ClausenFunctions.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ClausenFunctions.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ClausenFunctions_nomain.o src/ClausenFunctions.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ClausenFunctions.o ${OBJECTDIR}/src/ClausenFunctions_nomain.o;\
	fi

${OBJECTDIR}/src/BernoulliNumbers_nomain.o: ${OBJECTDIR}/src/BernoulliNumbers.o src/BernoulliNumbers.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/BernoulliNumbers.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I. -I. -I. -I. -I../gslpp/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/BernoulliNumbers_nomain.o src/BernoulliNumbers.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/BernoulliNumbers.o ${OBJECTDIR}/src/BernoulliNumbers_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f3 || true; \
	    ${TESTDIR}/TestFiles/f2 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libloopfunctions.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
