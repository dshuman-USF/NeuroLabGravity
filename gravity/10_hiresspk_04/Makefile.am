bin_PROGRAMS = hires

hires_SOURCES = hires_04.f colors32.f colors.f	\
plane_poly_04.f surfacev3.f			\
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

LDADD = $(x11_LIBS)  $(FLIBS) $(SB_LIBS) -lncurses # $(top_builddir)/common/libcommon.a

files: $(hires_SOURCES) Makefile.am
	ls $(hires_SOURCES) | grep -v common > files

AM_CFLAGS = --std=c99 -I$(top_builddir)/common $(SB_CFLAGS)
