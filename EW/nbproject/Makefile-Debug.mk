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
	${OBJECTDIR}/src/Alepton.o \
	${OBJECTDIR}/src/obliqueShat.o \
	${OBJECTDIR}/src/AFBcharm.o \
	${OBJECTDIR}/src/Acharm.o \
	${OBJECTDIR}/src/PtauPol.o \
	${OBJECTDIR}/src/Rlepton.o \
	${OBJECTDIR}/src/Mw.o \
	${OBJECTDIR}/src/obliqueS.o \
	${OBJECTDIR}/src/Rbottom.o \
	${OBJECTDIR}/src/obliqueY.o \
	${OBJECTDIR}/src/obliqueW.o \
	${OBJECTDIR}/src/obliqueT.o \
	${OBJECTDIR}/src/sigmaHadron.o \
	${OBJECTDIR}/src/obliqueThat.o \
	${OBJECTDIR}/src/GammaZ.o \
	${OBJECTDIR}/src/AFBbottom.o \
	${OBJECTDIR}/src/GammaW.o \
	${OBJECTDIR}/src/obliqueEpsilon1.o \
	${OBJECTDIR}/src/Rcharm.o \
	${OBJECTDIR}/src/sin2thetaEff.o \
	${OBJECTDIR}/src/obliqueEpsilon3.o \
	${OBJECTDIR}/src/EW.o \
	${OBJECTDIR}/src/obliqueU.o \
	${OBJECTDIR}/src/Abottom.o \
	${OBJECTDIR}/src/AFBlepton.o \
	${OBJECTDIR}/src/obliqueEpsilon2.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a

${OBJECTDIR}/src/Alepton.o: src/Alepton.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Alepton.o src/Alepton.cpp

${OBJECTDIR}/src/obliqueShat.o: src/obliqueShat.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueShat.o src/obliqueShat.cpp

${OBJECTDIR}/src/AFBcharm.o: src/AFBcharm.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AFBcharm.o src/AFBcharm.cpp

${OBJECTDIR}/src/Acharm.o: src/Acharm.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Acharm.o src/Acharm.cpp

${OBJECTDIR}/src/PtauPol.o: src/PtauPol.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/PtauPol.o src/PtauPol.cpp

${OBJECTDIR}/src/Rlepton.o: src/Rlepton.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Rlepton.o src/Rlepton.cpp

${OBJECTDIR}/src/Mw.o: src/Mw.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Mw.o src/Mw.cpp

${OBJECTDIR}/src/obliqueS.o: src/obliqueS.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueS.o src/obliqueS.cpp

${OBJECTDIR}/src/Rbottom.o: src/Rbottom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Rbottom.o src/Rbottom.cpp

${OBJECTDIR}/src/obliqueY.o: src/obliqueY.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueY.o src/obliqueY.cpp

${OBJECTDIR}/src/obliqueW.o: src/obliqueW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueW.o src/obliqueW.cpp

${OBJECTDIR}/src/obliqueT.o: src/obliqueT.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueT.o src/obliqueT.cpp

${OBJECTDIR}/src/sigmaHadron.o: src/sigmaHadron.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/sigmaHadron.o src/sigmaHadron.cpp

${OBJECTDIR}/src/obliqueThat.o: src/obliqueThat.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueThat.o src/obliqueThat.cpp

${OBJECTDIR}/src/GammaZ.o: src/GammaZ.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/GammaZ.o src/GammaZ.cpp

${OBJECTDIR}/src/AFBbottom.o: src/AFBbottom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AFBbottom.o src/AFBbottom.cpp

${OBJECTDIR}/src/GammaW.o: src/GammaW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/GammaW.o src/GammaW.cpp

${OBJECTDIR}/src/obliqueEpsilon1.o: src/obliqueEpsilon1.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueEpsilon1.o src/obliqueEpsilon1.cpp

${OBJECTDIR}/src/Rcharm.o: src/Rcharm.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Rcharm.o src/Rcharm.cpp

${OBJECTDIR}/src/sin2thetaEff.o: src/sin2thetaEff.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/sin2thetaEff.o src/sin2thetaEff.cpp

${OBJECTDIR}/src/obliqueEpsilon3.o: src/obliqueEpsilon3.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueEpsilon3.o src/obliqueEpsilon3.cpp

${OBJECTDIR}/src/EW.o: src/EW.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/EW.o src/EW.cpp

${OBJECTDIR}/src/obliqueU.o: src/obliqueU.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueU.o src/obliqueU.cpp

${OBJECTDIR}/src/Abottom.o: src/Abottom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/Abottom.o src/Abottom.cpp

${OBJECTDIR}/src/AFBlepton.o: src/AFBlepton.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/AFBlepton.o src/AFBlepton.cpp

${OBJECTDIR}/src/obliqueEpsilon2.o: src/obliqueEpsilon2.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -g -I../gslpp/src -I../Observables/src -I../StandardModel/src -I../EWSM/src -I../ZFitter/src -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/obliqueEpsilon2.o src/obliqueEpsilon2.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libew.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
