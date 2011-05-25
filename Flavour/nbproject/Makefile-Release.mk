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
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/Vcb.o \
	${OBJECTDIR}/_ext/1086379272/HeffDF2.o \
	${OBJECTDIR}/src/gamma.o \
	${OBJECTDIR}/src/Vus.o \
	${OBJECTDIR}/src/Vub.o \
	${OBJECTDIR}/src/alpha_2a.o \
	${OBJECTDIR}/src/alpha.o \
	${OBJECTDIR}/src/DmBd.o \
	${OBJECTDIR}/src/DmBd0.o \
	${OBJECTDIR}/src/Vud.o \
	${OBJECTDIR}/src/EvolDF2.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a

${OBJECTDIR}/src/Vcb.o: src/Vcb.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vcb.o src/Vcb.cpp

${OBJECTDIR}/_ext/1086379272/HeffDF2.o: /afs/infn.it/roma1/project/susy/susy/SusyFit/Flavour/src/HeffDF2.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1086379272
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/1086379272/HeffDF2.o /afs/infn.it/roma1/project/susy/susy/SusyFit/Flavour/src/HeffDF2.cpp

${OBJECTDIR}/src/gamma.o: src/gamma.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gamma.o src/gamma.cpp

${OBJECTDIR}/src/Vus.o: src/Vus.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vus.o src/Vus.cpp

${OBJECTDIR}/src/Vub.o: src/Vub.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vub.o src/Vub.cpp

${OBJECTDIR}/src/alpha_2a.o: src/alpha_2a.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/alpha_2a.o src/alpha_2a.cpp

${OBJECTDIR}/src/alpha.o: src/alpha.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/alpha.o src/alpha.cpp

${OBJECTDIR}/src/DmBd.o: src/DmBd.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DmBd.o src/DmBd.cpp

${OBJECTDIR}/src/DmBd0.o: src/DmBd0.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/DmBd0.o src/DmBd0.cpp

${OBJECTDIR}/src/Vud.o: src/Vud.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Vud.o src/Vud.cpp

${OBJECTDIR}/src/EvolDF2.o: src/EvolDF2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EvolDF2.o src/EvolDF2.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libflavour.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
