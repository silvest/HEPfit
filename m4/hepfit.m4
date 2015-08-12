dnl ###### GSL CHECK ######
AC_DEFUN([HEPFIT_CHECK_GSL],
[
AC_MSG_CHECKING([for gsl location])
GSLINCLUDE=""
GSLLIBS=""

AC_ARG_WITH(gsl, AC_HELP_STRING([],[]),
        [],
	[with_gsl=system])

if test "x$with_gsl" = "xno"; then
AC_MSG_ERROR([libgsl is required. Please install the GNU scientific library and header files.])
fi

if test "x$with_gsl" = "xsystem"; then
	AC_MSG_RESULT([in system libraries or non-standard directory])
	oldlibs="$LIBS"
	AC_CHECK_LIB(m,main)
	AC_CHECK_LIB(gslcblas,main)
	AC_CHECK_LIB(gsl,main,[with_gsl=check],
			[
			AC_MSG_ERROR([Cannot find libgsl. Please install the GNU scientific library and header files or use --with-gsl=.])
			]
		     )
	GSLLIBS="$LIBS"
	LIBS=$oldlibs

	
dnl	if test "x$with_gsl" = "xcheck"; then
dnl	
dnl	fi

else
	if test "`uname -m`" = "x86_64" -a -e "$with_gsl/lib64/libgsl.a" -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib64 -R$with_gsl/lib64 -lgslcblas -lgsl"
		GSLINCLUDE="-I$with_gsl/include"
	elif test -e "$with_gsl/lib/libgsl.a" -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib -R$with_gsl/lib -lgslcblas -lgsl"
		GSLINCLUDE="-I$with_gsl/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_gsl/lib/libgsl.a or the headers in $with_gsl/include])
	fi
fi

dnl AC_SUBST(GSLINCLUDE)
dnl AC_SUBST(GSLLIBS)
])

dnl ###### LOOPTOOLS PATH ######
AC_DEFUN([HEPFIT_LOOPTOOLS_PATH],
[
AC_MSG_CHECKING([for looptools location])

AC_ARG_WITH(looptools,
        AC_HELP_STRING([--with-looptools=DIR],[location of looptools installation @<:@default=system libs@:>@]),
        [],
	[with_looptools=system])

if test "x$with_looptools" = "xno"; then
AC_MSG_ERROR([liblooptools is required. Please install LoopTools.])
fi

if test "x$with_looptools" = "xsystem"; then
	AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_CHECK_LIB(ooptools,main,[],
			[
			AC_MSG_ERROR([Cannot find liblooptools. Please install LoopTools or use --with-looptools=.])
			]
		     )
	LTLIBS="$LIBS"
	LIBS=$oldlibs
else
	if test "`uname -m`" = "x86_64" -a -e "$with_looptools/lib64/libooptools.a" -a -d "$with_looptools/include/"; then
		AC_MSG_RESULT([found in $with_looptools])
		LTLIBS="-L$with_looptools/lib64 -R$with_looptools/lib64 -looptools"
		LTINCLUDE="-I$with_looptools/include"
	elif test -e "$with_looptools/lib/libooptools.a" -a -d "$with_looptools/include/"; then
		AC_MSG_RESULT([found in $with_looptools])
		LTLIBS="-L$with_looptools/lib -R$with_looptools/lib -looptools"
		LTINCLUDE="-I$with_looptools/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_looptools/lib/libooptools.a or the headers in $with_looptools/include])
	fi
fi

SuF_LIBS="$SuF_LIBS $LTLIBS"
SuF_CFLAGS="$SuF_CFLAGS $LTINCLUDE"

dnl # AC_SUBST(GSLINCLUDE)
dnl # AC_SUBST(GSLLIBS)
])


dnl ##### LOOPTOOLS #####
AC_DEFUN([HEPFIT_LOOPTOOLS],
[
dnl #AC_REQUIRE([AC_PROG_FC])
dnl #AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
dnl AC_REQUIRE([AC_PROG_F77])
dnl AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
dnl AC_REQUIRE([AC_PROG_CC([gcc])])


AC_MSG_CHECKING([if Looptools build works])
enable_looptools=yes

if test "x$GCC" = "xyes"; then
   case "${host}" in
      x86_64-*|*-darwin1*)
	AM_FCFLAGS="$AM_FCFLAGS -fdefault-integer-8"
      	;;
   esac

dnl    AC_LANG_PUSH([Fortran])
   	oldFCFLAGS="$FCFLAGS"
   	FCFLAGS="$AM_FCFLAGS"
   	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      print *[,]"Hello"])],
		[],
		[AC_MSG_RESULT([no])
 		 AC_MSG_ERROR([needs gfortran on 64bit machines])]
	)
	FCFLAGS="$oldFCFLAGS"
dnl   AC_LANG_POP([Fortran])
fi
AC_MSG_RESULT([$enable_looptools])

dnl # AC_SUBST([F77],[$FC])
dnl # AC_SUBST([FFLAGS],[$FCFLAGS])
dnl # AC_SUBST([AM_FFLAGS],[$AM_FCFLAGS])
dnl # AC_SUBST([FLIBS],[$FCLIBS])
])


dnl ###### BAT PATH ######
AC_DEFUN([HEPFIT_BAT_PATH],
[
AC_MSG_CHECKING([for BAT location])

AC_ARG_WITH(bat,
        AC_HELP_STRING([--with-bat=DIR],[location of BAT installation @<:@default=system libs@:>@]),
        [],
	[with_bat=system])

if test "x$with_bat" = "xno"; then
AC_MSG_ERROR([libBAT is required. Please install BAT.])
fi

if test "x$with_bat" = "xsystem"; then
	AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_CHECK_LIB(BAT,main,[],
			[
			AC_MSG_ERROR([Cannot find libBAT. Please install BAT or use --with-bat=.])
			]
		     )
	BATLIBS="$LIBS"
	LIBS=$oldlibs
else
	if test "`uname -m`" = "x86_64" -a -e "$with_bat/lib64/libBAT.a" -a -d "$with_looptools/include/BAT"; then
		AC_MSG_RESULT([found in $with_bat])
		BATLIBS="-L$with_bat/lib64 -R$with_bat/lib64 -lBATmodels -lBAT"
		BATINCLUDE="-I$with_bat/include"
	elif test -e "$with_bat/lib/libBAT.a" -a -d "$with_bat/include/BAT"; then
		AC_MSG_RESULT([found in $with_bat])
		BATLIBS="-L$with_bat/lib -R$with_bat/lib -lBATmodels -lBAT"
		BATINCLUDE="-I$with_bat/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_bat/lib/libBAT.a or the headers in $with_bat/include])
	fi
fi

SuF_LIBS="$SuF_LIBS $BATLIBS"
SuF_CFLAGS="$SuF_CFLAGS $BATINCLUDE"
])

dnl ###### FeynHiggs PATH ######
dnl # TODO: impose required version
dnl # TODO: check if compiled with gfortran, else return error message
AC_DEFUN([HEPFIT_FEYNHIGGS_PATH],
[
AC_MSG_CHECKING([for feynhiggs location])

AC_ARG_WITH(feynhiggs,
        AC_HELP_STRING([--with-feynhiggs=DIR],[location of FeynHiggs installation @<:@default=system libs@:>@]),
        [],
	[with_feynhiggs=system])

if test "x$with_feynhiggs" = "xno"; then
AC_MSG_ERROR([libFH is required. Please install FeynHiggs.])
fi

if test "x$with_feynhiggs" = "xsystem"; then
	AC_MSG_RESULT([in system libraries or non-standard directory])
	oldlibs="$LIBS"
	AC_CHECK_LIB(FH,main,[],
			[
			AC_MSG_ERROR([Cannot find libFH. Please install FeynHiggs or use --with-feynhiggs=.])
			]
		     )
	FHLIBS="$LIBS"
	LIBS=$oldlibs
else
	if test "`uname -m`" = "x86_64" -a -e "$with_feynhiggs/lib64/libFH.a" -a -d "$with_feynhiggs/include/"; then
		AC_MSG_RESULT([found in $with_feynhiggs])
		FHLIBS="-L$with_feynhiggs/lib64 -R$with_feynhiggs/lib64 -lFH"
		FHINCLUDE="-I$with_feynhiggs/include"
	elif test -e "$with_feynhiggs/lib/libFH.a" -a -d "$with_bat/include/"; then
		AC_MSG_RESULT([found in $with_feynhiggs])
		FHLIBS="-L$with_feynhiggs/lib -R$with_feynhiggs/lib -lFH"
		FHINCLUDE="-I$with_feynhiggs/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_feynhiggs/lib/libFH.a or the headers in $with_feynhiggs/include])
	fi
fi

SuF_LIBS="$SuF_LIBS $FHLIBS"
SuF_CFLAGS="$SuF_CFLAGS $FHINCLUDE"
])



dnl ## HEPFIT overview
dnl # TODO: is python really required? by root? and ..?
dnl #*** debug mode:	$enable_debug
dnl #***
AC_DEFUN([HEPFIT_OVERVIEW],
[
FCSTRING=`$FC --version | head -1`
CXXSTRING=`$CXX --version | head -1`
CCSTRING=`$CC --version | head -1`
if test "x$PYTHON" != "x:"
then
   python_was_found="yes, using Python $PYTHON_VERSION"
else
   python_was_found="no, requires Python >= 2.6"
fi
cat << _HEPFIT_EOF_ > config.hepfit
*****************************************************
*** $PACKAGE_STRING configuration summary
*** Please include this information in bug reports!
***--------------------------------------------------
*** Prefix:		$prefix
***
*** BAT:		$with_bat
*** FeynHiggs:		$with_feynhiggs
***
*** GSL:		`$GSL_CONFIG --prefix`
*** root:   		${ROOTCFLAGS} (root-version: $ROOTVERSION)
*** boost:              ${BOOST_CPPFLAGS}
***
*** Host:		$host
*** CC:			$CCSTRING
*** CXX:		$CXXSTRING
*** FC:			$FCSTRING
***
*** CXXFLAGS:		$CXXFLAGS
*****************************************************
_HEPFIT_EOF_
])
