SHELL		= /bin/sh

CXX		= g++
LINKER		= g++

OPT		= -O
PDBINCDIR	= ../libpdb++
CFLAGS		= $(OPT) -I$(PDBINCDIR)

LFLAGS		=
PDBLIBDIR	= ../libpdb++
LIBRARIES	= -L$(PDBLIBDIR) -lpdb++ -lm

#
# MODIFY ITEMS BELOW AT YOUR OWN RISK
#

.SUFFIXES:	.cc
.cc.o:
	$(CXX) $(CFLAGS) -c $< -o $@

PROG	= ksdssp

HDRS	= Atom.h ListArray.h Residue.h Structure.h misc.h \
	  List.h Model.h SquareArray.h ksdssp.h
SRCS	= ksdssp.cc Model.cc Residue.cc Structure.cc misc.cc
OBJS	= ksdssp.o Model.o Residue.o Structure.o misc.o

$(PROG):	$(OBJS)
	$(LINKER) $(LFLAGS) $(OBJS) $(LIBRARIES) -o $@

ksdssp.csh: ksdssp.csh.sed
	sed s\|BINDIR\|${BINDIR}\| < ${.CURDIR}/ksdssp.csh.sed > ${.TARGET}
	chmod +x ${.TARGET}

clean:
	-rm -f $(OBJS)
	-rm -rf ii_files cxx_repository

distclean:	clean
	-rm -f $(PROG)

ksdssp.o:	ksdssp.cc ksdssp.h Model.h Residue.h Atom.h List.h \
		ListArray.h SquareArray.h Structure.h \
		${PDBINCDIR}/pdb++.h

Model.o:	Model.cc ksdssp.h Model.h Residue.h Atom.h List.h \
		ListArray.h SquareArray.h Structure.h misc.h \
		${PDBINCDIR}/pdb++.h

Residue.o:	Residue.cc ksdssp.h Residue.h Atom.h List.h misc.h

Structure.o:	Structure.cc Structure.h List.h

misc.o:		misc.cc ksdssp.h misc.h 
