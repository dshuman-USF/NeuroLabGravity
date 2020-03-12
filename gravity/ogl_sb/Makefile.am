noinst_PROGRAMS = test_ogl_sb

test_ogl_sb_SOURCES = test_ogl_sb.f ogl_sb.c windata.h util.h	\
view_camera.c text3d.c oglstroke.c


AM_FFLAGS =  $(EXTRA_FFLAGS) -I$(top_builddir)/common
AM_CFLAGS = --std=c99 $(SB_CFLAGS) -I$(top_builddir)/common

LDADD = $(x11_LIBS) -lm $(FLIBS) -lGL $(SB_LIBS)

files: $(test_ogl_sb_SOURCES) Makefile.am
	ls $(test_ogl_sb_SOURCES) > files



