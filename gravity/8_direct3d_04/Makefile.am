bin_PROGRAMS = direct3d_bl direct3d_bl_mp direct3d_bl_sig_01	\
direct3d_bl_sig_05 direct3d_bl_sub

common_sources = colors32.f head_04.def		\
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

direct3d_bl_SOURCES        = direct3d_bl_04.f       $(common_sources)
direct3d_bl_mp_SOURCES     = direct3d_bl_mp_04_2.f  $(common_sources)
direct3d_bl_sig_01_SOURCES = direct3d_bl_sig_01_2.f $(common_sources)
direct3d_bl_sig_05_SOURCES = direct3d_bl_sig_05_2.f $(common_sources)
direct3d_bl_sub_SOURCES    = direct3d_bl_sub_04_2.f $(common_sources)

all_sources = $(direct3d_bl_SOURCES) $(direct3d_bl_mp_SOURCES)	\
$(direct3d_bl_sig_01_SOURCES) $(direct3d_bl_sig_05_SOURCES)	\
$(direct3d_bl_sub_SOURCES)

AM_FFLAGS =  $(EXTRA_FFLAGS) -I$(top_builddir)/common

LDADD = $(x11_LIBS)  $(FLIBS) $(SB_LIBS) -lncurses # $(top_builddir)/common/libcommon.a

files: $(all_sources) Makefile.am
	ls $(all_sources) | grep -v common | sort -u > files

AM_CFLAGS = --std=c99 -I$(top_builddir)/common $(SB_CFLAGS)
