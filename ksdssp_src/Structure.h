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

#ifndef structure_h
#define structure_h

#include "List.h"

#define	B_PARA	1
#define	B_ANTI	2

class Sheet;

class Helix {
	int	from_, to_;
	int	type_;		// For pdb HELIX record
public:
		Helix(int from, int to)
			{ from_ = from; to_ = to; type_ = 0; }
	int	type(void) const { return type_; }
	int	from(void) const { return from_; }
	int	to(void) const { return to_; }
	void	setType(int t) { type_ = t; }
};

class Ladder {
	char	name_;
	int	type_;
	int	start_[2], end_[2];
	Ladder	*neighbor_[2];
	Sheet	*sheet_;
	int	isBulge_;
public:
		Ladder(int type, int s1, int e1, int s2, int e2);
	int	type(void) const { return type_; }
	int	start(int n) const { return start_[n]; }
	int	end(int n) const { return end_[n]; }
	char	name(void) const { return name_; }
	void	setName(char n);
	Sheet	*sheet(void) const { return sheet_; }
	void	setSheet(Sheet *s) { sheet_ = s; }
	Ladder	*neighbor(int n) const { return neighbor_[n]; }
	void	setNeighbor(int n, Ladder *l) { neighbor_[n] = l; }
	int	isBulge(void) const { return isBulge_; }
	void	setBulge(void) { isBulge_ = 1; }
	int	overlaps(const Ladder *l, int overlap[2]) const;
	Ladder	*otherNeighbor(Ladder *l) const;
	int	neighborCount(void) const;
public:
	static Ladder
		*mergeBulge(const Ladder *l1, const Ladder *l2);
};

class Sheet {
	char		name_;
	List<Ladder>	ladderList_;
public:
		Sheet(char name) : ladderList_() { name_ = name; }
	char	name(void) const { return name_; }
	void	addLadder(Ladder *l) { ladderList_.append(l); }
	const List<Ladder> &
		ladderList(void) const { return ladderList_; }
	Ladder	*firstLadder(void);
};

#endif
