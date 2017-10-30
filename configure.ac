#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
AC_INIT([RFMIX], [v2.03-r0], markoni.wright@gmail.com)
AC_CONFIG_SRCDIR([./rfmix.cpp])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdint.h stdlib.h string.h errno.h sys/time.h unistd.h])
AC_CHECK_FUNCS(dup2 gettimeofday)

AC_CHECK_LIB([pthread], [pthread_create], [], AC_MSG_FAILURE([POSIX thread library (-lpthread) is not available. Can not build or run without this.]),[-lpthread])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gettimeofday memset strdup strerror strtol strtoul])

AC_OUTPUT(Makefile)
