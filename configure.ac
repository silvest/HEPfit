AC_PREREQ([2.59])
AC_INIT([HEPfit], [1.0], [hepfit@roma1.infn.it], [HEPfit], [hepfit.roma1.infn.it])
AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([Config/config.h])

AM_INIT_AUTOMAKE
AM_PROG_AR

AC_LANG(C++)
AM_INIT_AUTOMAKE([1.10 -Wall -Werror no-define])
LT_INIT


AC_PROG_CC([gcc])
AC_PROG_CXX([g++])
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

AX_PATH_GSL(1.11, , AC_MSG_ERROR(could not find required version of GSL)) #GSL 
ROOT_PATH(3.03/05, AC_DEFUN([HAVE_ROOT]), AC_MSG_ERROR(ROOT is not found or version is too old)) #ROOT
HEPFIT_CHECK_GSL
HEPFIT_LOOPTOOLS_PATH
BOOST_REQUIRE
BOOST_PROGRAM_OPTIONS  
HEPFIT_BAT_PATH


AC_SUBST(SuF_CFLAGS)
AC_SUBST(SuF_LIBS)

AC_SUBST([HEPfit_SO_VERSION], [1:0:0])
AC_SUBST([HEPfit_API_VERSION], [1.0])

AC_CONFIG_FILES([Makefile
		HEPfit-${HEPfit_API_VERSION}.pc:HEPfit.pc.in])

AC_CONFIG_FILES([hepfit-config], [chmod +x hepfit-config])

AC_SUBST([hepfitCflags],["$(export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:.; ./hepfit-config --cflags)"])
AC_SUBST([hepfitLibs],["$(export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:.; ./hepfit-config --libs)"])

HEPFIT_OVERVIEW

AC_CONFIG_COMMANDS([summary],[cat config.hepfit])

AC_OUTPUT
