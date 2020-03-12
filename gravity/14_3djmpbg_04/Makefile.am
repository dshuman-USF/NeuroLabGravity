bin_PROGRAMS = 3djmp

3djmp_SOURCES = 3djmp.f colors32.f rana2.f	\
pace.c                                          \
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

AM_FFLAGS = $(EXTRA_FFLAGS) -I$(top_builddir)/common	\
-ffixed-line-length-0

LDADD = $(x11_LIBS)  $(FLIBS) $(SB_LIBS) -lncurses # $(top_builddir)/common/libcommon.a

files: $(3djmp_sources) Makefile.am
	ls $(3djmp_sources) | grep -v common | sort -u > files

AM_CFLAGS = --std=c99 -I$(top_builddir)/common -I$(top_builddir)/ogl_sb $(SB_CFLAGS)

3djmp.o: colors32.o
