/*
 * Copyright (c) 2002 The Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions, and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions, and the following
 *      disclaimer in the documentation and/or other materials provided
 *      with the distribution.
 *   3. Redistributions must acknowledge that this software was
 *      originally developed by the UCSF Computer Graphics Laboratory
 *      under support by the NIH National Center for Research Resources,
 *      grant P41-RR01081.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 * OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef residue_h
#define residue_h

#ifdef NEED_BOOL
#define	bool	int
#define	false	0
#define	true	1
#endif

#include <stdio.h>
#include <pdb++.h>
#include "Atom.h"
#include "List.h"

#define	R_3DONOR	0x0001
#define	R_3ACCEPTOR	0x0002
#define	R_3GAP		0x0004
#define	R_3HELIX	0x0008
#define	R_4DONOR	0x0010
#define	R_4ACCEPTOR	0x0020
#define	R_4GAP		0x0040
#define	R_4HELIX	0x0080
#define	R_PBRIDGE	0x0100
#define	R_ABRIDGE	0x0200
#define	R_TER		0x8000

class Residue;

class Residue {
	PDB::Residue	residue_;
	List<Atom>	aList_;
	int		flags_;
public:
			Residue(const PDB::Residue &r);
			~Residue(void);
	const PDB::Residue &
			residue(void) const { return residue_; }
	void		addAtom(Atom *a);
	int		addImideHydrogen(const Residue *prev);
	Atom		*atom(const std::string &name) const;
	int		sameAs(const PDB::Residue &r) const;
	int		hBondedTo(const Residue *other) const;
	int		printAtoms(FILE *output, int sn) const;
	void		printSummary(FILE *output) const;
	int		flag(int f) const;
	void		setFlag(int f);
public:
	static void	setHBondCutoff(float cutoff);
};

inline
Residue::Residue(const PDB::Residue &r)
	: aList_()
{
	residue_ = r;
	flags_ = 0;
}

inline
Residue::~Residue(void)
{
	for (Pix p = aList_.first(); p != 0; aList_.next(p))
		delete aList_(p);
}

inline void
Residue::addAtom(Atom *a)
{
	aList_.append(a);
}

inline Atom *
Residue::atom(const std::string &name) const
{
	for (Pix p = aList_.first(); p != 0; aList_.next(p)) {
		Atom *a = aList_(p);
		if (a->name() == name)
			return a;
	}
	return NULL;
}

inline int
Residue::flag(int f) const
{
	return flags_ & f;
}

inline void
Residue::setFlag(int f)
{
	flags_ |= f;
}

#endif
