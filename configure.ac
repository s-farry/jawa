AC_INIT([Jawa], [1.0], [bug-automake@gnu.org])

dnl AC_CONFIG_MACRO_DIR([/user2/sfarry/Jawa/m4])

m4_include([m4/root.m4])
dnl m4_include([m4/ax_boost_base.m4])
dnl m4_include([m4/ax_boost_python.m4])
dnl m4_include([m4/boost.m4])
dnl m4_include([m4/ax_python.m4])
m4_include([m4/ax_python_devel.m4])
m4_include([m4/ax_boost_base.m4])
m4_include([m4/ax_boost_python.m4])

AC_PROG_CXX
dnl AX_PYTHON
AX_PYTHON_DEVEL
dnl BOOST_REQUIRE
dnl AX_BOOST_PYTHON

ROOT_PATH(5.34/07, , AC_MSG_ERROR(Failed: Your ROOT version is too old))
ROOT_FEATURE([python], ,[AC_MSG_ERROR(Failed: No python with root install)])

AC_OUTPUT(Makefile)
