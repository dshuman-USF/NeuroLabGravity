bin_PROGRAMS = blkshuf

blkshuf_SOURCES = blkshuf_05.f  rdfrm5.f ransub.f

AM_FFLAGS =  $(EXTRA_FFLAGS)

LDADD = $(FLIBS)

files: $(blkshuf_SOURCES) Makefile.am
	ls $(blkshuf_SOURCES) > files
