bin_PROGRAMS = xgtune

xgtune_SOURCES = xgtune_04.f distunev4.f	\
$(top_builddir)/common/sbparam.defs		\
$(top_builddir)/common/sbparam.h		\
$(top_builddir)/common/xsbf.c			\
$(top_builddir)/common/ps_starbase.c		\
$(top_builddir)/common/ps_starbase.h		\
$(top_builddir)/common/config.defs

AM_FFLAGS =  $(EXTRA_FFLAGS) -I$(top_builddir)/common
AM_CFLAGS =  $(EXTRA_CFLAGS) -Wall --std=c99

LDADD = $(x11_LIBS) -lm $(FLIBS) $(top_builddir)/common/libcommon.a

files: $(xgtune_SOURCES) Makefile.am
	ls $(xgtune_SOURCES) > files
