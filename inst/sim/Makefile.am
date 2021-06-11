lib_LTLIBRARIES = libmcsim.la

libmcsim_la_SOURCES = \
	delays.c \
	delays.h \
	getopt.c \
	getopt.h \
	hungtype.h \
	lex.c lex.h \
	lexerr.c \
	lexerr.h \
	lexfn.c \
	lexfn.h \
	list.c \
	list.h \
	lsodes.h \
	lsodes1.c \
	lsodes2.c \
	matutil.c \
	matutil.h \
	matutilo.c \
	matutilo.h \
	mh.c \
	mh.h \
	modelu.c \
	modelu.h \
	modiface.h \
	optdsign.c \
	optdsign.h \
	random.c \
	random.h \
	sim.c \
	sim.h \
	simi.c \
	simi.h \
	siminit.c \
	siminit.h \
	simmonte.c \
	simmonte.h \
	simo.c \
	simo.h \
	strutil.c \
	strutil.h \
	yourcode.c \
	yourcode.h

include_HEADERS = \
	delays.h \
	matutil.h \
	list.h \
	lex.h \
	sim.h \
	lexfn.h \
	hungtype.h \
	modiface.h \
	modelu.h \
	random.h \
	yourcode.h

AM_CFLAGS = -O3 -w

simdir=$(datadir)/mcsim/sim
sim_DATA = $(libmcsim_la_SOURCES)

bin_SCRIPTS = makemcsim makemcsims makemcsimd makemcsimp

CLEANFILES = $(bin_SCRIPTS)
