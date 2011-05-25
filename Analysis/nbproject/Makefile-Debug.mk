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
CND_PLATFORM=GNU-Linux-x86
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
LDLIBSOPTIONS=-L/usr/lib/root -L../Flavour/dist/Debug/GNU-Linux-x86 -L../gslpp/dist/Debug/GNU-Linux-x86 -L../Observables/dist/Debug/GNU-Linux-x86 -L../InputParser/dist/Debug/GNU-Linux-x86 -L../StandardModel/dist/Debug/GNU-Linux-x86 -L../Utils/dist/Debug/GNU-Linux-x86 -L../MonteCarlo/dist/Debug/GNU-Linux-x86 -L../BAT/lib -lgsl -lgslcblas -lBAT -lBATmodels -lboost_program_options ../MonteCarlo/dist/Debug/GNU-Linux-x86/libmontecarlo.a ../InputParser/dist/Debug/GNU-Linux-x86/libinputparser.a ../Flavour/dist/Debug/GNU-Linux-x86/libflavour.a ../Observables/dist/Debug/GNU-Linux-x86/libobservables.a ../StandardModel/dist/Debug/GNU-Linux-x86/libstandardmodel.a ../Utils/dist/Debug/GNU-Linux-x86/libutils.a ../gslpp/dist/Debug/GNU-Linux-x86/libgslpp.a -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lMinuit -lThread -lGui -lm -ldl

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../MonteCarlo/dist/Debug/GNU-Linux-x86/libmontecarlo.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../InputParser/dist/Debug/GNU-Linux-x86/libinputparser.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Flavour/dist/Debug/GNU-Linux-x86/libflavour.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Observables/dist/Debug/GNU-Linux-x86/libobservables.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../StandardModel/dist/Debug/GNU-Linux-x86/libstandardmodel.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../Utils/dist/Debug/GNU-Linux-x86/libutils.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ../gslpp/dist/Debug/GNU-Linux-x86/libgslpp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -pthread -rdynamic -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../Observables/src -I../MonteCarlo/src -I../gslpp/src -I../Flavour/src -I../StandardModel/src -I../InputParser/src -I../BAT/include -I/usr/include/ -I/usr/include/root -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/analysis

# Subprojects
.clean-subprojects:
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Flavour && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../gslpp && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../InputParser && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../MonteCarlo && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Observables && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../StandardModel && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
