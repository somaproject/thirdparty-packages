dnl Copyright (c) 2009 by CodeSourcery, Inc.  All rights reserved.
dnl
dnl File:   release.m4
dnl Author: Stefan Seefeld
dnl Date:   2009-03-26
dnl
dnl Contents: release-related configuration for Sourcery VSIPL++
dnl

AC_DEFUN([SVXX_RELEASE],
[

# Activate internal CodeSourcery documentation makefiles.
AC_ARG_ENABLE(csl-documentation,,
  [case "$enableval" in
    yes) enable_csl_documentation=1 ;;
    *)  enable_csl_documentation= ;;
   esac])
AC_SUBST(enable_csl_documentation)

# Specify the csl-docbook directory, used for internal documentation.
AC_ARG_WITH(csl-docbook,,
  [csl_docbook_prefix=$withval],
  [csl_docbook_prefix="$srcdir/doc/csl-docbook"])

if test "$enable_csl_documentation" = 1; then
  AC_CHECK_FILE($csl_docbook_prefix/GNUmakefile.inc,,
    [AC_MSG_ERROR([Please set csl-docbook prefix via --with-csl-docbook.])])
fi
AC_SUBST(csl_docbook_prefix)

# Specify the xml catalog files to use when building internal documentation.
AC_ARG_WITH(xml-catalog-files,,
  [case "$withval" in
     yes) AC_MSG_ERROR([xml catalog files not specified]) ;;
     no) ;;
     *)   XML_CATALOG_FILES="$withval" ;;
   esac])
AC_SUBST(XML_CATALOG_FILES)

AC_ARG_WITH(version-string,
  AS_HELP_STRING([--with-version-string=VERSION],
		 [The version string, with just the version number.]),
  [case "$withval" in
    yes) AC_MSG_ERROR([version-string not specified]) ;;
    no) ;;
    *) version_string="$withval" 
       major_version_string="${withval%-*}"
       ;;
  esac])
if test "$enable_csl_documentation" = 1; then
  test "$version_string" || AC_MSG_ERROR([--with-version-string is required])
fi
AC_SUBST(version_string)
AC_SUBST(major_version_string)

AC_ARG_WITH(pkgversion,
  AS_HELP_STRING([--with-pkgversion=VERSION],
		 [The version string, including package prefixes.]),
  [case "$withval" in
    yes) AC_MSG_ERROR([pkgversion not specified]) ;;
    no) ;;
    *) pkgversion="$withval" ;;
  esac])
if test "$enable_csl_documentation" = 1; then
  test "$pkgversion" || AC_MSG_ERROR([--with-pkgversion is required])
fi
AC_SUBST(pkgversion)

])
