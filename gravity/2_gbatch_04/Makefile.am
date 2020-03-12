bin_PROGRAMS = gbatch

gbatch_SOURCES = gbatch.f gmint.f possub.f icel_win.f postodir3d.f	\
ghead.def $(top_builddir)/common/config.defs

AM_FFLAGS = $(EXTRA_FFLAGS) -I$(top_builddir)/common	\
-ffixed-line-length-0
AM_CFLAGS =  $(EXTRA_CFLAGS) -Wall --std=c99

LDADD = $(x11_LIBS) -lm $(FLIBS)

files: $(gbatch_SOURCES) Makefile.am
	ls $(gbatch_SOURCES) > files
