
dnl Initialize the build system
dnl DONT FORGET TO CHANGE THE VERSION!
AC_INIT([CCwcs],[0.0.2],[jvinniec@gmail.edu])
AC_CONFIG_SRCDIR([CCwcs/])

dnl Run the automake stuff
#AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror])

# Make sure there is a c++ compiler installed
AC_PROG_CC
AC_PROG_CXX([clang++ g++])
AC_LANG_PUSH([C++])

dnl Other required programs
dnl AC_PROG_RANLIB
AM_PROG_AR
LT_INIT

dnl This variable will be used specify special compile options
special_exec=

CCWCS_CPPFLAGS=

dnl #################################
dnl ADDITIONAL OPTIONS FOR '--enable-'/'--disable-'
dnl #################################


dnl #################################
dnl ADDITIONAL OPTIONS FOR '--with-'/'--without-'
dnl #################################



dnl #################################
dnl CHECK FOR SPECIFIED PROGRAMS
dnl #################################


dnl #################################
dnl CHECK FOR STANDARD STUFF
dnl #################################

dnl ---------------
dnl option checks
dnl ---------------


dnl ---------------
dnl library checks
dnl ---------------

dnl Check for standard library functions
#AC_CHECK_FUNCS([floor mkdir pow sqrt])
#AC_C_INLINE
#AC_TYPE_SIZE_T
#AC_TYPE_UINT32_T

dnl ---------------
dnl header checks
dnl ---------------

dnl Some initial checks
AC_CHECK_HEADER_STDBOOL

dnl C++ headers
#AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])


dnl #################################
dnl SETUP WCS
dnl #################################

dnl ---------------
dnl option checks
dnl ---------------

dnl WCS options
wcslibs=
WCS_CXXFLAGS=
have_wcs=no

AC_ARG_WITH([wcs],
	AS_HELP_STRING([--with-wcs],
		[directory containing 'libwcs/libwcs.a' and 'libwcs/wcs.h']))
AC_ARG_WITH([wcs-libdir],
	AS_HELP_STRING([--with-wcs-libdir],
		[directory containing 'libwcs.a']))
AC_ARG_WITH([wcs-incdir],
	AS_HELP_STRING([--with-wcs-incdir=PATH],
		[directory containing 'wcs.h']))

dnl The following is true if '--with-wcs' is specified
if (test "x$with_wcs" != x) && (test -d ${with_wcs}/libwcs); then
   wcslibs="-L$with_wcs/libwcs ${with_wcs}/libwcs/libwcs.a -lwcs $wcslibs"
   WCS_CXXFLAGS="-I$with_wcs/libwcs ${WCS_CXXFLAGS}"
fi

dnl The following is true if '--with-wcs-include' is provided
if test "x$with_wcs_incdir" != x; then
   WCS_CXXFLAGS="-I$with_wcs_incdir ${WCS_CXXFLAGS}"
fi

dnl The following is true if '--with-wcs-lib' is provided
if test "x$with_wcs_libdir" != x; then
   wcslibs="-L$with_wcs_libdir $wcslibs"
fi

dnl The following is true if none of the '--with-wcs...' options
dnl were specified
#if test "x$with_wcs$with_wcs_include$with_wcs_lib" = x ; then
#fi

dnl Add the include flags and the library search paths
LIBS="${wcslibs} ${LIBS}"


dnl ---------------
dnl library checks
dnl ---------------

have_wcs=no
AC_CHECK_LIB([wcs],[wcsdist],
	[have_wcs=yes],
	[AC_MSG_ERROR([wcs library does not contain 'wcsdist'])])
AC_CHECK_LIB([wcs],[wcskinit],
	[have_wcs=yes],
	[AC_MSG_ERROR([wcs library does not contain 'wcskinit'])])
AC_CHECK_LIB([wcs],[wcscon],
	[have_wcs=yes],
	[AC_MSG_ERROR([wcs library does not contain 'wcscon'])])
AC_CHECK_LIB([wcs],[tnxpix],
	[have_wcs=yes],
	[AC_MSG_ERROR([wcs library does not contain 'tnxpix'])])
AC_CHECK_LIB([wcs],[tnxpos],
	[have_wcs=yes],
	[AC_MSG_ERROR([wcs library does not contain 'tnxpos'])])
dnl If we have a valid version of wcs, add the flag to the library
if test "x${have_wcs}" != xyes; then
   echo "------------------------------------"
   echo " ERROR :: -lwcs"
   echo " Could not find required methods in"
   echo " wcs library."
   echo "------------------------------------"
   (exit 1); exit 1;
fi
if test "x${have_wcs}" = xyes; then
   LIBS="-lwcs $LIBS"
fi

dnl ---------------
dnl header checks
dnl ---------------

dnl Check for wcs.h
dnl AC_CHECK_HEADERS([wcs.h])
AC_SUBST(WCS_CXXFLAGS)


dnl #################################
dnl ADD THE CCWCS_CPPFLAGS
dnl #################################
AC_SUBST(CCWCS_CPPFLAGS)

dnl
AC_CONFIG_HEADERS([config.h])

dnl #################################
dnl Specify the Makefile locations
dnl #################################
AC_CONFIG_FILES([Makefile CCwcs/Makefile])

dnl
AC_OUTPUT

dnl Some helpful information to the user
echo \
"
=================================================
 ${PACKAGE_NAME} Version ${PACKAGE_VERSION}
=================================================
 Prefix: '${prefix}'
 C Compiler   : ${CC} 
 C++ Compiler : ${CXX}
 CFLAGS  : ${CFLAGS} 
 CPPFLAGS: ${CPPFLAGS}
 CXXFLAGS: ${CXXFLAGS}

 Additional compiler flags:
 --------------------------
${WCS_CXXFLAGS} ${CCWCS_CPPFLAGS}

 LDFLAGS : ${LDFLAGS}

 Additional linker flags:
 ------------------------
${wcslibs} $LIBS

 Package features:
   Special Execution: ${special_exec}

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all             - build all binaries
     install         - install everything
=================================================
"