AC_DEFUN([ACX_SCALAPACK],
[
acx_scalapack_ok=no

AC_ARG_WITH(scalapack,
	[AC_HELP_STRING([--with-scalapack=<lib>], [use ScaLAPACK library <lib>])])
case $with_scalapack in
	yes | "") ;;
	no) acx_scalapack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) SCALAPACK_LIBS="$with_scalapack" ;;
	*) SCALAPACK_LIBS="-l$with_scalapack" ;;
esac

acx_scalapack_save_LIBS="$LIBS"

# First, check SCALAPACK_LIBS environment variable
if test "x$SCALAPACK_LIBS" != x; then
# Get fortran symbol name of pdpotrf (a ScaLAPACK routine).
	pdpotrf=pdpotrf_
	save_LIBS="$LIBS"; LIBS="$SCALAPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $pdpotrf in $SCALAPACK_LIBS])
	AC_TRY_LINK_FUNC($pdpotrf, [acx_scalapack_ok=yes], [SCALAPACK_LIBS=""])
	AC_MSG_RESULT($acx_scalapack_ok)
	LIBS="$save_LIBS"
    AC_DEFINE(HAVE_SCALAPACK, "1")
fi
# restore the libs variable
LIBS=$acx_scalapack_save_LIBS
]) 
