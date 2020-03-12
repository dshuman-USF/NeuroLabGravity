bin_PROGRAMS = spkpat6bg spkpat6kbg spkpat6bgr spkpatwip spkpatdist svddemo

spkpat6bg_SOURCES = surdat.f spkpat6bg_04.f   smooth1.f pcg_basic.c pcg_basic.h
spkpat6kbg_SOURCES = surdat.f spkpat6kbg_04.f smooth1.f pcg_basic.c pcg_basic.h
spkpat6bgr_SOURCES = surdat.f spkpat6bgr.f    smooth1.f pcg_basic.c pcg_basic.h ranz.f
spkpatwip_SOURCES = spkpat.f smooth1.f
spkpatdist_SOURCES = spkpatdist.f pcg_basic.c pcg_basic.h
svddemo_SOURCES = svddemo.f pcg_basic.c pcg_basic.h

spkpat6bg.o : surdat.o

spkpat6kbg.o : surdat.o

AM_FFLAGS = $(EXTRA_FFLAGS) -ffixed-line-length-0 $(LAPACK_CFLAGS)

LDADD = $(FLIBS) $(LAPACK_LIBS)

files: $(spkpat6bg_SOURCES) $(spkpat6kbg_SOURCES) Makefile.am
	ls $(spkpat6bg_SOURCES) $(spkpat6kbg_SOURCES) | sort -u > files
