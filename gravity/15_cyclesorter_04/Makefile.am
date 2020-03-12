bin_PROGRAMS = cyclesorter

cyclesorter_SOURCES = cyclesorter_04.f pick_04.f sort.f

AM_FFLAGS =  $(EXTRA_FFLAGS)

LDADD = $(FLIBS)

files: $(cyclesorter_SOURCES) Makefile.am
	ls $(cyclesorter_SOURCES) > files
