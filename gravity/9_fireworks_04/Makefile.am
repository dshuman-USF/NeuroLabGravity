bin_PROGRAMS = fireworks

fireworks_SOURCES = fireworks_04.f colors32.f rana2.f index.f	\
$(top_builddir)/common/sbparam.defs				\
$(top_builddir)/common/sbparam.h				\
$(top_builddir)/common/xsbf.c					\
$(top_builddir)/common/ps_starbase.c				\
$(top_builddir)/common/ps_starbase.h				\
$(top_builddir)/common/config.defs

AM_FFLAGS =  $(EXTRA_FFLAGS) -I$(top_builddir)/common

LDADD = $(x11_LIBS) $(top_builddir)/common/libcommon.a

files: $(fireworks_SOURCES) Makefile.am
	ls $(fireworks_SOURCES) | grep -v common > files

AM_CFLAGS = --std=c99
