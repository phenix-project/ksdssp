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

#include <math.h>
#include <string.h>
#include "ksdssp.h"
#include "Residue.h"
#include "misc.h"

#ifndef DONT_INSTANIATE
template class List<Atom>;
#endif

static float	hBondCutoff = -0.5;

//
// Check if two PDB residues are the same
//
int
Residue::sameAs(const PDB::Residue &r) const
{
	if (r.seqNum != residue_.seqNum)
		return 0;
	if (r.chainId != residue_.chainId)
		return 0;
	if (r.insertCode != residue_.insertCode)
		return 0;
	if (strncmp(r.name, residue_.name, sizeof (PDB::RName)) != 0)
		return 0;
	return 1;
}

//
// Add the imide hydrogen if it is missing
//
int
Residue::addImideHydrogen(const Residue *prev)
{
	if (prev == NULL || atom(" H") != NULL)
		return 0;		// Already there
	Atom *n = atom(" N");
	if (n == NULL) {
		if (verbose)
			(void) fprintf(stderr,
				"N missing in residue %d%c[%c]\n",
				residue_.seqNum, residue_.chainId,
				residue_.insertCode);
		return -1;
	}
	Atom *ca = atom(" CA");
	if (ca == NULL) {
		if (verbose)
			(void) fprintf(stderr,
				"CA missing in residue %d%c[%c]\n",
				residue_.seqNum, residue_.chainId,
				residue_.insertCode);
		return -1;
	}
	Atom *c = prev->atom(" C");
	if (c == NULL) {
		if (verbose)
			(void) fprintf(stderr,
				"C missing in residue %d%c[%c]\n",
				prev->residue().seqNum,
				prev->residue().chainId,
				prev->residue().insertCode);
		return -1;
	}
	Atom *o = prev->atom(" O");
	if (o == NULL) {
		if (verbose)
			(void) fprintf(stderr,
				"O missing in residue %d%c[%c]\n",
				prev->residue().seqNum,
				prev->residue().chainId,
				prev->residue().insertCode);
		return -1;
	}

	const float *nCoord = n->coord();
	const float *caCoord = ca->coord();
	const float *cCoord = c->coord();
	const float *oCoord = o->coord();
	float v1[3], v2[3], v3[3];
	int i;
	for (i = 0; i < 3; i++) {
		v1[i] = caCoord[i] - nCoord[i];
		v2[i] = cCoord[i] - nCoord[i];
		v3[i] = oCoord[i] - cCoord[i];
	}
	normalize(v1);
	normalize(v2);
	normalize(v3);
	float p1[3], hDir[3];
	bisect(p1, v1, v2);
	bisect(hDir, p1, v3);

	const float nhLength = 1.01;
	float hCoord[3];
	for (i = 0; i < 3; i++)
		hCoord[i] = nCoord[i] - nhLength * hDir[i];

	addAtom(new Atom(" H", hCoord));
	return 0;
}

//
// Check if other residue is hydrogen bonded to this one
//
int
Residue::hBondedTo(const Residue *other) const
{
	const float q1 = 0.42;
	const float q2 = 0.20;
	const float f = 332;

	Atom *c = atom(" C");
	Atom *o = atom(" O");
	if (c == NULL || o == NULL)
		return 0;
	Atom *n = other->atom(" N");
	Atom *h = other->atom(" H");
	if (n == NULL || h == NULL)
		return 0;
	float rCN = distSquared(c->coord(), n->coord());
	if (rCN > 49.0)		// Optimize a little bit
		return 0;
	rCN = sqrtf(rCN);
	float rON = distance(o->coord(), n->coord());
	float rCH = distance(c->coord(), h->coord());
	float rOH = distance(o->coord(), h->coord());

	float E = q1 * q2 * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN) * f;
	return E < hBondCutoff;
}

//
// Print atom list to output stream
//
int
Residue::printAtoms(FILE *output, int sn) const
{
	for (Pix p = aList_.first(); p != 0; aList_.next(p)) {
		Atom *a = aList_(p);
		const float *c = a->coord();
		PDB pdb(PDB::ATOM);
		PDB::Atom &atom = pdb.atom;
		atom.serialNum = sn++;
		(void) strncpy(atom.name, a->name().c_str(),
				sizeof (PDB::AName));
		atom.residue = residue_;
		for (int i = 0; i < 3; i++)
			atom.xyz[i] = c[i];
		(void) fprintf(output, "%s\n", pdb.chars());
	}
	return sn;
}

//
// Print summary of residue state
//
void
Residue::printSummary(FILE *output) const
{
	char summary = ' ';
	if (flag(R_3HELIX))
		summary = 'G';
	else if (flag(R_4HELIX))
		summary = 'H';
	else if (flag(R_PBRIDGE | R_ABRIDGE))
		summary = 'E';

	char turn3 = ' ';
	if (flag(R_3DONOR) && flag(R_3ACCEPTOR))
		turn3 = 'X';
	else if (flag(R_3ACCEPTOR))
		turn3 = '>';
	else if (flag(R_3DONOR))
		turn3 = '<';
	else if (flag(R_3GAP))
		turn3 = '3';

	char turn4 = ' ';
	if (flag(R_4DONOR) && flag(R_4ACCEPTOR))
		turn4 = 'X';
	else if (flag(R_4ACCEPTOR))
		turn4 = '>';
	else if (flag(R_4DONOR))
		turn4 = '<';
	else if (flag(R_4GAP))
		turn4 = '4';

	char bridge = ' ';
	if (flag(R_PBRIDGE) && flag(R_ABRIDGE))
		bridge = '+';
	else if (flag(R_PBRIDGE))
		bridge = 'p';
	else if (flag(R_ABRIDGE))
		bridge = 'A';

	(void) fprintf(output, "%4.4s %4d%c[%c] -> %c %c %c %c\n",
		residue_.name, residue_.seqNum,
		residue_.chainId, residue_.insertCode,
		summary, turn3, turn4, bridge);
}

//
// Set the energy cutoff for hydrogen bondedness
//
void
Residue::setHBondCutoff(float cutoff)
{
	hBondCutoff = cutoff;
}
