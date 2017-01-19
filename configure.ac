AC_INIT([Jawa], [1.0], [bug-automake@gnu.org])

## Check and block installation into the src/build dir
if test "$prefix" = "$PWD"; then
  AC_MSG_ERROR([Installation into the build directory is not supported: use a different --prefix argument])
fi
  ## Force default prefix to have a path value rather than NONE
if test "$prefix" = "NONE"; then
   prefix=/usr/local
fi

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CXX
LT_INIT([disable-static])
AM_DISABLE_STATIC
AM_PROG_LIBTOOL
AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE
AC_ARG_ENABLE([python],
[  --enable-python    Enable Python],
[case "${enableval}" in
  yes) python=true ;;
  no)  python=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-python]) ;;
esac],[debug=false])
AM_CONDITIONAL([WANT_PYTHON], [test x$python = xtrue])


## Extend and substitute the default build flags after lib testing
AM_CPPFLAGS="-I\$(top_srcdir)/include -I\$(top_builddir)/include"
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_CXXFLAGS)

AC_PROG_CXX

ROOT_PATH(6.00/00, , AC_MSG_ERROR(Failed: Your ROOT version is too old))
if test x$python = xtrue; then
   dnl BOOST_REQUIRE
   AX_PYTHON_DEVEL	 
   AX_BOOST_PYTHON
   ROOT_FEATURE([python], ,[AC_MSG_ERROR(Failed: No python with root install)])
fi

AC_OUTPUT(Makefile src/Makefile include/Makefile jawa.pc)

