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
	${OBJECTDIR}/src/TwoLoopQCD.o \
	${OBJECTDIR}/src/TwoLoopEW.o \
	${OBJECTDIR}/src/ThreeLoopEW.o \
	${OBJECTDIR}/src/EWSM.o \
	${OBJECTDIR}/src/OneLoopEW.o \
	${OBJECTDIR}/src/EWSMcommon.o \
	${OBJECTDIR}/src/ApproximateFormulae.o \
	${OBJECTDIR}/src/ThreeLoopEW2QCD.o \
	${OBJECTDIR}/src/ThreeLoopQCD.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a

${OBJECTDIR}/src/TwoLoopQCD.o: src/TwoLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/TwoLoopQCD.o src/TwoLoopQCD.cpp

${OBJECTDIR}/src/TwoLoopEW.o: src/TwoLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/TwoLoopEW.o src/TwoLoopEW.cpp

${OBJECTDIR}/src/ThreeLoopEW.o: src/ThreeLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopEW.o src/ThreeLoopEW.cpp

${OBJECTDIR}/src/EWSM.o: src/EWSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSM.o src/EWSM.cpp

${OBJECTDIR}/src/OneLoopEW.o: src/OneLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/OneLoopEW.o src/OneLoopEW.cpp

${OBJECTDIR}/src/EWSMcommon.o: src/EWSMcommon.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMcommon.o src/EWSMcommon.cpp

${OBJECTDIR}/src/ApproximateFormulae.o: src/ApproximateFormulae.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ApproximateFormulae.o src/ApproximateFormulae.cpp

${OBJECTDIR}/src/ThreeLoopEW2QCD.o: src/ThreeLoopEW2QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopEW2QCD.o src/ThreeLoopEW2QCD.cpp

${OBJECTDIR}/src/ThreeLoopQCD.o: src/ThreeLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopQCD.o src/ThreeLoopQCD.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/testclass.o ${TESTDIR}/tests/testrunner.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}  -L/usr/lib/root -L/usr/local/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} -lcppunit -lgsl -lgslcblas ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a ../StandardModel/dist/Debug/GNU-MacOSX/libstandardmodel.a ../LoopFunctions/dist/Debug/GNU-MacOSX/libloopfunctions.a 


${TESTDIR}/tests/testclass.o: tests/testclass.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -Isrc -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/testclass.o tests/testclass.cpp


${TESTDIR}/tests/testrunner.o: tests/testrunner.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -Isrc -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -MMD -MP -MF $@.d -o ${TESTDIR}/tests/testrunner.o tests/testrunner.cpp


${OBJECTDIR}/src/TwoLoopQCD_nomain.o: ${OBJECTDIR}/src/TwoLoopQCD.o src/TwoLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/TwoLoopQCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/TwoLoopQCD_nomain.o src/TwoLoopQCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/TwoLoopQCD.o ${OBJECTDIR}/src/TwoLoopQCD_nomain.o;\
	fi

${OBJECTDIR}/src/TwoLoopEW_nomain.o: ${OBJECTDIR}/src/TwoLoopEW.o src/TwoLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/TwoLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/TwoLoopEW_nomain.o src/TwoLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/TwoLoopEW.o ${OBJECTDIR}/src/TwoLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/ThreeLoopEW_nomain.o: ${OBJECTDIR}/src/ThreeLoopEW.o src/ThreeLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ThreeLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopEW_nomain.o src/ThreeLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ThreeLoopEW.o ${OBJECTDIR}/src/ThreeLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/EWSM_nomain.o: ${OBJECTDIR}/src/EWSM.o src/EWSM.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSM.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSM_nomain.o src/EWSM.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSM.o ${OBJECTDIR}/src/EWSM_nomain.o;\
	fi

${OBJECTDIR}/src/OneLoopEW_nomain.o: ${OBJECTDIR}/src/OneLoopEW.o src/OneLoopEW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/OneLoopEW.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/OneLoopEW_nomain.o src/OneLoopEW.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/OneLoopEW.o ${OBJECTDIR}/src/OneLoopEW_nomain.o;\
	fi

${OBJECTDIR}/src/EWSMcommon_nomain.o: ${OBJECTDIR}/src/EWSMcommon.o src/EWSMcommon.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/EWSMcommon.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EWSMcommon_nomain.o src/EWSMcommon.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/EWSMcommon.o ${OBJECTDIR}/src/EWSMcommon_nomain.o;\
	fi

${OBJECTDIR}/src/ApproximateFormulae_nomain.o: ${OBJECTDIR}/src/ApproximateFormulae.o src/ApproximateFormulae.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ApproximateFormulae.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ApproximateFormulae_nomain.o src/ApproximateFormulae.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ApproximateFormulae.o ${OBJECTDIR}/src/ApproximateFormulae_nomain.o;\
	fi

${OBJECTDIR}/src/ThreeLoopEW2QCD_nomain.o: ${OBJECTDIR}/src/ThreeLoopEW2QCD.o src/ThreeLoopEW2QCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ThreeLoopEW2QCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopEW2QCD_nomain.o src/ThreeLoopEW2QCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ThreeLoopEW2QCD.o ${OBJECTDIR}/src/ThreeLoopEW2QCD_nomain.o;\
	fi

${OBJECTDIR}/src/ThreeLoopQCD_nomain.o: ${OBJECTDIR}/src/ThreeLoopQCD.o src/ThreeLoopQCD.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/ThreeLoopQCD.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I../gslpp/src -I../StandardModel/src -I../LoopFunctions/src -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThreeLoopQCD_nomain.o src/ThreeLoopQCD.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/ThreeLoopQCD.o ${OBJECTDIR}/src/ThreeLoopQCD_nomain.o;\
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
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
