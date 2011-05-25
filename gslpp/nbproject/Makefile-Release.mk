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
	${OBJECTDIR}/src/gslpp_vector_double.o \
	${OBJECTDIR}/src/gslpp_matrix_complex.o \
	${OBJECTDIR}/src/gslpp_complex.o \
	${OBJECTDIR}/src/gslpp_vector_complex.o \
	${OBJECTDIR}/src/gslpp_matrix_double.o

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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a

${OBJECTDIR}/src/gslpp_vector_double.o: src/gslpp_vector_double.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_vector_double.o src/gslpp_vector_double.cpp

${OBJECTDIR}/src/gslpp_matrix_complex.o: src/gslpp_matrix_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_matrix_complex.o src/gslpp_matrix_complex.cpp

${OBJECTDIR}/src/gslpp_complex.o: src/gslpp_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_complex.o src/gslpp_complex.cpp

${OBJECTDIR}/src/gslpp_vector_complex.o: src/gslpp_vector_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_vector_complex.o src/gslpp_vector_complex.cpp

${OBJECTDIR}/src/gslpp_matrix_double.o: src/gslpp_matrix_double.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} $@.d
	$(COMPILE.cc) -O2 -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_matrix_double.o src/gslpp_matrix_double.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/gslpptest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/f2: ${TESTDIR}/tests/newsimpletest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f2 $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/gslpptest.o: tests/gslpptest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${TESTDIR}/tests/gslpptest.o tests/gslpptest.cpp


${TESTDIR}/tests/newsimpletest.o: tests/newsimpletest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} $@.d
	$(COMPILE.cc) -O2 -I. -I. -I. -I. -MMD -MP -MF $@.d -o ${TESTDIR}/tests/newsimpletest.o tests/newsimpletest.cpp


${OBJECTDIR}/src/gslpp_vector_double_nomain.o: ${OBJECTDIR}/src/gslpp_vector_double.o src/gslpp_vector_double.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gslpp_vector_double.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_vector_double_nomain.o src/gslpp_vector_double.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gslpp_vector_double.o ${OBJECTDIR}/src/gslpp_vector_double_nomain.o;\
	fi

${OBJECTDIR}/src/gslpp_matrix_complex_nomain.o: ${OBJECTDIR}/src/gslpp_matrix_complex.o src/gslpp_matrix_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gslpp_matrix_complex.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_matrix_complex_nomain.o src/gslpp_matrix_complex.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gslpp_matrix_complex.o ${OBJECTDIR}/src/gslpp_matrix_complex_nomain.o;\
	fi

${OBJECTDIR}/src/gslpp_complex_nomain.o: ${OBJECTDIR}/src/gslpp_complex.o src/gslpp_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gslpp_complex.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_complex_nomain.o src/gslpp_complex.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gslpp_complex.o ${OBJECTDIR}/src/gslpp_complex_nomain.o;\
	fi

${OBJECTDIR}/src/gslpp_vector_complex_nomain.o: ${OBJECTDIR}/src/gslpp_vector_complex.o src/gslpp_vector_complex.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gslpp_vector_complex.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_vector_complex_nomain.o src/gslpp_vector_complex.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gslpp_vector_complex.o ${OBJECTDIR}/src/gslpp_vector_complex_nomain.o;\
	fi

${OBJECTDIR}/src/gslpp_matrix_double_nomain.o: ${OBJECTDIR}/src/gslpp_matrix_double.o src/gslpp_matrix_double.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gslpp_matrix_double.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} $@.d;\
	    $(COMPILE.cc) -O2 -Dmain=__nomain -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/gslpp_matrix_double_nomain.o src/gslpp_matrix_double.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gslpp_matrix_double.o ${OBJECTDIR}/src/gslpp_matrix_double_nomain.o;\
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
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgslpp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
