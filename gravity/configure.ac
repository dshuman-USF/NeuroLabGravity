#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.56)

AC_INIT(gravity, 0.17.11, dshuman@.usf.edu, gravity)

AM_INIT_AUTOMAKE(subdir-objects)

AC_CONFIG_SRCDIR([configure.ac])

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_F77(f90 gfortran)
if test x$F77 = xgfortran; then
   AC_SUBST(EXTRA_FFLAGS,["-O3 -fno-underscoring -fno-backslash -Wall -fcray-pointer -frecord-marker=4 -fconvert=big-endian -fno-automatic"])
fi
AC_PROG_LN_S

# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS
AC_CANONICAL_BUILD
case $build in
     *hp*) 
           AC_SUBST(x11_LIBS,["-L/opt/graphics/common -L/opt/graphics/common/lib -lXwindow -lhpgfx -lXhp11 -lX11"])
           AC_SUBST(EXTRA_FFLAGS,["-K +U77"])
           AC_SUBST(EXTRA_CFLAGS,[-mgas]);;
     *)    PKG_CHECK_MODULES([x11],[x11]);;
esac
PKG_CHECK_MODULES(GSL, [gsl])
PKG_CHECK_MODULES(SB, [ftgl libpng glu])
PKG_CHECK_MODULES(LAPACK, [lapack])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF(void *)
SIZEOF_VOID_P=$ac_cv_sizeof_void_p
AC_SUBST(SIZEOF_VOID_P)

# Checks for library functions.
AC_CONFIG_FILES([Makefile common/Makefile common/config.defs])
AC_CONFIG_FILES(            [2_gbatch_04/Makefile])
AC_CONFIG_FILES(           [4_xtrydis_04/Makefile])
AC_CONFIG_FILES(           [5_xprojtm_04/Makefile])
AC_CONFIG_FILES(            [6_xslope_09/Makefile])
AC_CONFIG_FILES(          [8_direct3d_04/Makefile])
AC_CONFIG_FILES(         [9_fireworks_04/Makefile])
AC_CONFIG_FILES(         [10_hiresspk_04/Makefile])
AC_CONFIG_FILES(            [11_gtune_04/Makefile])
AC_CONFIG_FILES(     [12a_random_shuf_05/Makefile])
AC_CONFIG_FILES([12_local_block_shift_04/Makefile])
AC_CONFIG_FILES(           [13_spkpat_04/Makefile])
AC_CONFIG_FILES(          [14_3djmpbg_04/Makefile])
AC_CONFIG_FILES(      [15_cyclesorter_04/Makefile])
AC_CONFIG_FILES(                 [ogl_sb/Makefile])

AC_OUTPUT
