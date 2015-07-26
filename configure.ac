AC_INIT([Jawa], [1.0], [bug-automake@gnu.org])

dnl AC_CONFIG_MACRO_DIR([/user2/sfarry/Jawa/m4])

m4_include([m4/root.m4])
m4_include([m4/ax_boost_base.m4])
m4_include([m4/ax_boost_python.m4])
m4_include([m4/ax_python.m4])
m4_include([m4/ax_python_devel.m4])

AC_PROG_CXX
AX_PYTHON
AX_PYTHON_DEVEL

ROOT_PATH(5.34/07, , AC_MSG_ERROR(Your ROOT version is too old))
ROOT_FEATURE([python], ,[AC_MSG_ERROR(No python with root)])

AC_OUTPUT(Makefile)
