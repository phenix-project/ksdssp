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

#include <stdio.h>
#include <ctype.h>
#include "Structure.h"

//
// Constructor for Ladder
//
Ladder::Ladder(int type, int s1, int e1, int s2, int e2)
{
	name_ = '?';
	type_ = type;
	neighbor_[0] = NULL;
	neighbor_[1] = NULL;
	sheet_ = NULL;
	isBulge_ = 0;
	if (s1 < e1) {
		start_[0] = s1;
		end_[0] = e1;
	}
	else {
		start_[0] = e1;
		end_[0] = s1;
	}
	if (s2 < e2) {
		start_[1] = s2;
		end_[1] = e2;
	}
	else {
		start_[1] = e2;
		end_[1] = s2;
	}
}

//
// Set the name of this ladder
//
void
Ladder::setName(char n)
{
	if (type() == B_PARA)
		name_ = isupper(n) ? tolower(n) : n;
	else
		name_ = islower(n) ? toupper(n) : n;
}

//
// Check if two ladders overlap (share a residue)
//
int
Ladder::overlaps(const Ladder *l, int overlap[2]) const
{
	for (int i = 0; i < 2; i++) {
		int s1 = start(i);
		int e1 = end(i);
		for (int j = 0; j < 2; j++) {
			int s2 = l->start(j);
			int e2 = l->end(j);
			if (e1 >= s2 && e2 >= s1) {
				overlap[0] = i;
				overlap[1] = j;
				return 1;
			}
		}
	}
	return 0;
}

//
// Return neighbor which is not the given neighbor
//
Ladder *
Ladder::otherNeighbor(Ladder *l) const
{
	if (neighbor_[0] == l)
		return neighbor_[1];
	return neighbor_[0];
}

//
// Check number of neighbors
//
int
Ladder::neighborCount(void) const
{
	int count = 0;
	if (neighbor_[0] != NULL)
		count++;
	if (neighbor_[1] != NULL)
		count++;
	return count;
}

//
// Check whether two ladders should merge to form a beta bulge
// We take advantage of some properties of how the ladders were generated:
//	start/end(0) < start/end(1)
//
// Beta-bulge as defined by K&S:
//	"a bulge-linked ladder consists of two (perfect)
//	ladders or bridges of the same type connected by
//	at most one extra residue on one strand and at most
//	four residues on the other strand."
//
Ladder *
Ladder::mergeBulge(const Ladder *l1, const Ladder *l2)
{
	if (l1->type() != l2->type())
		return NULL;
	// Make sure that l1 precedes l2
	if (l1->start(0) > l2->start(0)) {
		const Ladder *tmp = l1;
		l1 = l2;
		l2 = tmp;
	}

	int d0 = l2->start(0) - l1->end(0);
	if (d0 < 0 || d0 > 4)
		return NULL;
	int d1;
	if (l1->type() == B_PARA)
		d1 = l2->start(1) - l1->end(1);
	else
		d1 = l1->start(1) - l2->end(1);
	if (d1 < 0 || d1 > 4)
		return NULL;
	if (d0 > 1 && d1 > 1)
		return NULL;

	int s0 = l1->start(0);
	int e0 = l2->end(0);
	int s1, e1;
	if (l1->type() == B_PARA) {
		s1 = l1->start(1);
		e1 = l2->end(1);
	}
	else {
		s1 = l2->start(1);
		e1 = l1->end(1);
	}
	Ladder *l = new Ladder(l1->type(), s0, e0, s1, e1);
	l->setBulge();
	return l;
}

//
// Find a ladder on one end of the sheet (or a random one
// if its a barrel)
//
Ladder *
Sheet::firstLadder(void)
{
	for (Pix p = ladderList_.first(); p != 0; ladderList_.next(p)) {
		Ladder *l = ladderList_(p);
		if (l->neighborCount() == 1)
			return l;
	}
	return ladderList_.head();
}
