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
	${OBJECTDIR}/src/ThreeLoopEW2QCD.o \
	${OBJECTDIR}/src/ThreeLoopQCD.o


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

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libewsm.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
