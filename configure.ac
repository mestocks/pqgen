AC_PREREQ([2.69])
AC_INIT([pqgen], [0.3.2], [mestocks@protonmail.com])

AC_PROG_CC
AC_PROG_INSTALL
AC_TYPE_UINT64_T
AC_FUNC_MALLOC

LT_INIT

AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([src/main.c])

AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h])

AC_CHECK_LIB([m],[exp],,[])

AM_INIT_AUTOMAKE([subdir-objects])

AC_OUTPUT




