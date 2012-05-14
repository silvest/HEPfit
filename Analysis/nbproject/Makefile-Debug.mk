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
	${OBJECTDIR}/main.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

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
LDLIBSOPTIONS=-L/usr/lib/root -L../Flavour/dist/Debug/GNU-Linux-x86 -L../EW/dist/Debug/GNU-Linux-x86 -L../gslpp/dist/Debug/GNU-Linux-x86 -L../Observables/dist/Debug/GNU-Linux-x86 -L../SUSYMassInsertion/dist/Debug/GNU-Linux-x86 -L../InputParser/dist/Debug/GNU-Linux-x86 -L../StandardModel/dist/Debug/GNU-Linux-x86 -L../Utils/dist/Debug/GNU-Linux-x86 -L../MonteCarlo/dist/Debug/GNU-Linux-x86 -L../BAT/lib -L../LoopFunctions/dist/Debug/GNU-Linux-x86 -L../MFV/dist/Debug/GNU-Linux-x86 -L../SUSY/dist/Debug/GNU-Linux-x86 -L../FeynHiggs/lib64 -lgsl -lgslcblas -lboost_program_options ../MonteCarlo/dist/Debug/GNU-MacOSX/libmontecarlo.a ../InputParser/dist/Debug/GNU-MacOSX/libinputparser.a ../Flavour/dist/Debug/GNU-MacOSX/libflavour.a ../EW/dist/Debug/GNU-MacOSX/libew.a ../LoopFunctions/dist/Debug/GNU-MacOSX/libloopfunctions.a ../Observables/dist/Debug/GNU-MacOSX/libobservables.a ../MFV/dist/Debug/GNU-MacOSX/libmfv.a ../SUSY/dist/Debug/GNU-MacOSX/libsusy.a ../StandardModel/dist/Debug/GNU-MacOSX/libstandardmodel.a ../Utils/dist/Debug/GNU-MacOSX/libutils.a ../SUSYMassInsertion/dist/Debug/GNU-MacOSX/libsusymassinsertion.a ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMinuit -lThread -lGui -lm -ldl -lgfortran -lBAT -lBATmodels -lFH

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../MonteCarlo/dist/Debug/GNU-MacOSX/libmontecarlo.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../InputParser/dist/Debug/GNU-MacOSX/libinputparser.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Flavour/dist/Debug/GNU-MacOSX/libflavour.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../EW/dist/Debug/GNU-MacOSX/libew.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../LoopFunctions/dist/Debug/GNU-MacOSX/libloopfunctions.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Observables/dist/Debug/GNU-MacOSX/libobservables.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../MFV/dist/Debug/GNU-MacOSX/libmfv.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../SUSY/dist/Debug/GNU-MacOSX/libsusy.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../StandardModel/dist/Debug/GNU-MacOSX/libstandardmodel.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Utils/dist/Debug/GNU-MacOSX/libutils.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../SUSYMassInsertion/dist/Debug/GNU-MacOSX/libsusymassinsertion.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../gslpp/dist/Debug/GNU-MacOSX/libgslpp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -pthread -rdynamic -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I/usr/include/ -I/usr/include/root -I../Utils/src -I../Observables/src -I../MonteCarlo/src -I../gslpp/src -I../Flavour/src -I../StandardModel/src -I../InputParser/src -I../BAT/include -I../EW/src -I../SUSYMassInsertion/src -I../LoopFunctions/src -I../MFV/src -I../SUSY/src -I../FeynHiggs/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug
	cd ../EW && ${MAKE}  -f Makefile CONF=Debug
	cd ../LoopFunctions && ${MAKE}  -f Makefile CONF=Debug
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug
	cd ../MFV && ${MAKE}  -f Makefile CONF=Debug
	cd ../SUSY && ${MAKE}  -f Makefile CONF=Debug
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug
	cd ../SUSYMassInsertion && ${MAKE}  -f Makefile CONF=Debug
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/RunningMass.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/RunningMass.o: tests/RunningMass.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -g -I. -I. -I. -I. -I/usr/include/ -I/usr/include/root -I../Utils/src -I../Observables/src -I../MonteCarlo/src -I../gslpp/src -I../Flavour/src -I../StandardModel/src -I../InputParser/src -I../BAT/include -I../EW/src -I../SUSYMassInsertion/src -I../LoopFunctions/src -I../MFV/src -I../SUSY/src -I../FeynHiggs/include -MMD -MP -MF $@.d -o ${TESTDIR}/tests/RunningMass.o tests/RunningMass.cpp


${OBJECTDIR}/main_nomain.o: ${OBJECTDIR}/main.o main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/main.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -g -I/usr/include/ -I/usr/include/root -I../Utils/src -I../Observables/src -I../MonteCarlo/src -I../gslpp/src -I../Flavour/src -I../StandardModel/src -I../InputParser/src -I../BAT/include -I../EW/src -I../SUSYMassInsertion/src -I../LoopFunctions/src -I../MFV/src -I../SUSY/src -I../FeynHiggs/include -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/main_nomain.o main.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/main.o ${OBJECTDIR}/main_nomain.o;\
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
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis

# Subprojects
.clean-subprojects:
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../EW && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../LoopFunctions && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../MFV && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../SUSY && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../SUSYMassInsertion && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
