bin_PROGRAMS = xprojtm

xprojtm_SOURCES = xprojtm_04.f colors32.f	\
head_04.def					\
$(top_builddir)/common/sbparam.defs		\
$(top_builddir)/common/sbparam.h		\
$(top_builddir)/ogl_sb/ogl_sb.c			\
$(top_builddir)/ogl_sb/text3d.c			\
$(top_builddir)/ogl_sb/oglstroke.c		\
$(top_builddir)/ogl_sb/view_camera.c		\
$(top_builddir)/ogl_sb/ogl_sb.h			\
$(top_builddir)/ogl_sb/windata.h		\
$(top_builddir)/ogl_sb/util.h			\
$(top_builddir)/common/config.defs		\
$(top_builddir)/common/mode.c

AM_FFLAGS =  $(EXTRA_FFLAGS) -I$(top_builddir)/common

LDADD = $(x11_LIBS) -lm $(FLIBS) $(SB_LIBS) -lncurses #$(top_builddir)/common/libcommon.a

files: $(xprojtm_SOURCES) Makefile.am
	ls $(xprojtm_SOURCES) | grep -v common > files

AM_CFLAGS = --std=c99 -I$(top_builddir)/common $(SB_CFLAGS)
