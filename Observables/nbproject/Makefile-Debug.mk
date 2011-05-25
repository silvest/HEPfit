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
	${OBJECTDIR}/src/Observable2D.o \
	${OBJECTDIR}/src/Observable.o \
	${OBJECTDIR}/src/ModelParameter.o \
	${OBJECTDIR}/src/Likelihood.o \
	${OBJECTDIR}/src/GaussianLikelihood.o \
	${OBJECTDIR}/src/ThObservable.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a

${OBJECTDIR}/src/Observable2D.o: src/Observable2D.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Observable2D.o src/Observable2D.cpp

${OBJECTDIR}/src/Observable.o: src/Observable.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Observable.o src/Observable.cpp

${OBJECTDIR}/src/ModelParameter.o: src/ModelParameter.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ModelParameter.o src/ModelParameter.cpp

${OBJECTDIR}/src/Likelihood.o: src/Likelihood.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Likelihood.o src/Likelihood.cpp

${OBJECTDIR}/src/GaussianLikelihood.o: src/GaussianLikelihood.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/GaussianLikelihood.o src/GaussianLikelihood.cpp

${OBJECTDIR}/src/ThObservable.o: src/ThObservable.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../Utils/src -I../gslpp/src -I../StandardModel/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/ThObservable.o src/ThObservable.cpp

# Subprojects
.build-subprojects:
	cd /afs/infn.it/roma1/project/susy/susy/SusyFit/gslpp && ${MAKE}  -f Makefile CONF=Debug
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug
	cd /afs/infn.it/roma1/project/susy/susy/SusyFit/StandardModel && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libobservables.a

# Subprojects
.clean-subprojects:
	cd /afs/infn.it/roma1/project/susy/susy/SusyFit/gslpp && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../Utils && ${MAKE}  -f Makefile CONF=Debug clean
	cd /afs/infn.it/roma1/project/susy/susy/SusyFit/StandardModel && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
