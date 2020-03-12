bin_PROGRAMS = xslope surstats xslopewip

noinst_PROGRAMS = expb

xslope_SOURCES = xslope_09_v31.f colors.f modej.f		\
plot_op2_v2.f plot_op3_v2.f indexx.f index2.f head_04.def	\
$(top_builddir)/common/sbparam.defs				\
$(top_builddir)/common/sbparam.h				\
$(top_builddir)/common/xsbf.c					\
$(top_builddir)/common/ps_starbase.c				\
$(top_builddir)/common/ps_starbase.h				\
$(top_builddir)/common/config.defs

xslopewip_SOURCES = xslope.f colors.f modej.f			\
plot_op2_v2.f plot_op3_v2.f indexx.f index2.f head_04.def	\
$(top_builddir)/common/sbparam.defs				\
$(top_builddir)/common/sbparam.h				\
$(top_builddir)/common/xsbf.c					\
$(top_builddir)/common/ps_starbase.c				\
$(top_builddir)/common/ps_starbase.h				\
$(top_builddir)/common/config.defs

surstats_SOURCES: surstats.c
expb_SOURCES: expb.f

AM_FFLAGS = $(EXTRA_FFLAGS) -I$(top_builddir)/common	\
-ffixed-line-length-0

AM_CFLAGS = --std=c99 -Wall $(GSL_CFLAGS)

$(xslope_OBJECTS): head_04.def

LDADD = $(x11_LIBS) -lm $(FLIBS) $(top_builddir)/common/libcommon.a $(GSL_LIBS)

files: $(xslope_SOURCES) Makefile.am
	ls $(xslope_SOURCES) | grep -v common > files
