AC_INIT(rhd_dq_module, version-0.1)
AC_ARG_VAR(ARCH, ARCH)
if test "$ARCH" = "" -a "$COBOLD_3DPARTY" = "no"
then
	AC_MSG_ERROR([ARCH variable must be defined.])
fi
if test "$ARCH" = "" -o "$COBOLD_3DPARTY" = "yes"
then
	DEPSRCPATH=3dparty/cobold
	F90_BASEPATH=$DEPSRCPATH
	COBOLD_3DPARTY="yes"
else
	DEPSRCPATH="$F90_BASEPATH/hd/rhd/$ARCH"
	COBOLD_3DPARTY="no"
	grep -q "function xkaros" $DEPSRCPATH/opta_routines.f90
	if test $? -ne 0
	then
		XKAROS_INTERFACE=-DXKAROS_INTERFACE
	fi
fi
AC_ARG_WITH([pyvisit],
	    AS_HELP_STRING([--with-pyvisit],
			   [Compile visualization tools for VisIt]),
			   AS_IF([test "x$VISIT_PATH" = "x"],
				 AC_MSG_ERROR([VISIT_PATH variable must be defined.])))
AS_IF([test "x$VISIT_PATH" != "x"],[VISIT_LDFLAG=-L$(dirname `find $VISIT_PATH/. -name "libsimV2.a"`)])
AC_ARG_VAR(COBOLD_3DPARTY, Use CO5BOLD 3d-party sources)
AC_ARG_VAR(F90_BASEPATH, CO5BOLD path)
AC_ARG_VAR(VISIT_PATH, Path to VisIt)
AC_MSG_CHECKING([CO5BOLD sources in $DEPSRCPATH])
if test -f $DEPSRCPATH/.configurerc
then
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
	AC_MSG_ERROR([Cannot find CO5BOLD '.configurerc' file. Check $F90_BASEPATH is properly set.])
fi
if [[[ "$(grep 'F90_PARALLEL=' $DEPSRCPATH/.configurerc | sed 's/F90_PARALLEL=//g')" =~ ^\"?(O|o)(P|p)(E|e)(N|n)(M|m)(P|p)\"?$ ]]]
then
	OPENMP="yes"
else
	OPENMP="no"
fi
FC=$(grep 'F90_COMPILER' $DEPSRCPATH/.configurerc | sed 's/F90_COMPILER=//g' | sed 's/^\"//g' | sed 's/\"$//g')
AC_PROG_FC()
AC_PROG_CC()
AC_OPENMP()
AC_SUBST(FC,$FC)
AC_SUBST(G77,$G77)
AC_SUBST(ARCH,$ARCH)
AC_SUBST(DEPSRCPATH, $DEPSRCPATH)
AC_SUBST(OPENMP, $OPENMP)
AC_SUBST(XKAROS_INTERFACE,$XKAROS_INTERFACE)
AC_SUBST(VISIT_LDFLAG,$VISIT_LDFLAG)
AC_SUBST(DEBUG,$DEBUG)
AC_OUTPUT(Makefile)
