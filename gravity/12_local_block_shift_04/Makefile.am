bin_PROGRAMS = lbs

lbs_SOURCES = lbs_04_v3.f  rdfrm5.f 

AM_FFLAGS =  $(EXTRA_FFLAGS)

LDADD = $(FLIBS)

files: $(lbs_SOURCES) Makefile.am
	ls $(lbs_SOURCES) > files
