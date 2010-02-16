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

#include <ctype.h>
#include <string.h>
#include "ksdssp.h"
#include "Model.h"
#include "misc.h"

#ifndef DONT_INSTANIATE
template class List<Residue>;
template class ListArray<Residue>;
template class SquareArray<char>;
template class List<Helix>;
template class List<Ladder>;
template class List<Sheet>;
#endif

static int	curModelNumber = -1;
static int	minStrandLength = 3;
static int	minHelixLength = 3;
static int	checkBulges = 1;

inline int
min(int i1, int i2)
{
	return i1 < i2 ? i1 : i2;
}

inline int
max(int i1, int i2)
{
	return i1 > i2 ? i1 : i2;
}

//
// Constructor for Model (read residues/atoms from PDB file)
//
Model::Model(FILE *input)
	: error_(), rList_(), helixList_(), ladderList_(), sheetList_(),
	  fileRecord_(PDB::USER_FILE)
{
	Residue *r = NULL;
	anyMore_ = 0;
	char buf[256];
	modelNumber_ = curModelNumber;
	while (fgets(buf, sizeof buf, input) != NULL) {
		PDB pdb(buf);
		switch (pdb.type()) {
		  case PDB::ATOM: {
			PDB::Atom &a = pdb.atom;
			if (a.residue.chainId == '\0')
				a.residue.chainId = ' ';
			if (a.residue.insertCode == '\0')
				a.residue.insertCode = ' ';
			if (r == NULL || !r->sameAs(a.residue)) {
				r = new Residue(a.residue);
				rList_.append(r);
			}
			r->addAtom(new Atom(a.name, a.xyz));
			break;
		  }
		  case PDB::TER:
			if (r != NULL)
				r->setFlag(R_TER);
			anyMore_ = 1;
			goto done;
		  case PDB::USER_FILE:
			fileRecord_ = pdb;
			curModelNumber = pdb.userFile.model;
			modelNumber_ = curModelNumber;
			break;
		  case PDB::END:
			curModelNumber = -1;
			anyMore_ = 1;
			goto done;
		}
	}
done:
	if (rList_.count() > 0) {
		rArray_ = new ListArray<Residue>(rList_);
		hBond_ = new SquareArray<char>(rArray_->count());
	}
}

//
// Destructor for Model (get rid of residue list)
//
Model::~Model(void)
{
	Pix p;
	for (p = rList_.first(); p != 0; rList_.next(p))
		delete rList_(p);
	for (p = helixList_.first(); p != 0; helixList_.next(p))
		delete helixList_(p);
	for (p = ladderList_.first(); p != 0; ladderList_.next(p))
		delete ladderList_(p);
	for (p = sheetList_.first(); p != 0; sheetList_.next(p))
		delete sheetList_(p);
	if (rList_.count() > 0) {
		delete rArray_;
		delete hBond_;
	}
}

//
// Define secondary structure using K&S definition
//
void
Model::defineSecondaryStructure(void)
{
	addImideHydrogens();
	findHBonds();

	findTurns(3);
	markHelices(3);
	findTurns(4);
	markHelices(4);
	findHelices();

	findBridges();
	findSheets();
}

//
// Print list of residues to given file
//
void
Model::printResidues(FILE *output) const
{
	int sn = 1;
	for (Pix p = rList_.first(); p != 0; rList_.next(p))
		sn = rList_(p)->printAtoms(output, sn);
}

//
// Print summary of residues
//
void
Model::printSummary(FILE *output) const
{
	(void) fputs("Helix Summary\n", output);
	Pix p;
	for (p = helixList_.first(); p != 0; helixList_.next(p)) {
		Helix *h = helixList_(p);
		const PDB::Residue &from = residue(h->from())->residue();
		const PDB::Residue &to = residue(h->to())->residue();
		(void) fprintf(output, "%2d: %4d%c[%c] -> %4d%c[%c]\n",
			h->type(),
			from.seqNum, from.chainId, from.insertCode,
			to.seqNum, to.chainId, to.insertCode);
	}
	(void) fputs("\n", output);

	(void) fputs("Ladder Summary\n", output);
	for (p = ladderList_.first(); p != 0; ladderList_.next(p)) {
		Ladder *l = ladderList_(p);
		const PDB::Residue &f0 = residue(l->start(0))->residue();
		const PDB::Residue &t0 = residue(l->end(0))->residue();
		const PDB::Residue &f1 = residue(l->start(1))->residue();
		const PDB::Residue &t1 = residue(l->end(1))->residue();
		(void) fprintf(output, "%c %4d%c[%c] -> %4d%c[%c] "
			"%-12s %4d%c[%c] -> %4d%c[%c]\n",
			l->name(),
			f0.seqNum, f0.chainId, f0.insertCode,
			t0.seqNum, t0.chainId, t0.insertCode,
			l->type() == B_PARA ? "parallel" : "antiparallel",
			f1.seqNum, f1.chainId, f1.insertCode,
			t1.seqNum, t1.chainId, t1.insertCode);
	}
	(void) fputs("\n", output);

	(void) fputs("Sheet Summary\n", output);
	for (p = sheetList_.first(); p != 0; sheetList_.next(p)) {
		Sheet *s = sheetList_(p);
		(void) fprintf(output, "Sheet %c:\n", s->name());
		Ladder *fl = s->firstLadder();
		Ladder *pl = NULL;
		Ladder *spl = NULL;
		for (Ladder *l = fl; l != NULL && !(l == fl && pl != NULL);
		spl = l, l = l->otherNeighbor(pl), pl = spl) {
			char n0 = l->neighbor(0) == NULL ? '-' :
					l->neighbor(0)->name();
			char n1 = l->neighbor(1) == NULL ? '-' :
					l->neighbor(1)->name();
			(void) fprintf(output, "\tLadder %c: %c %c\n",
						l->name(), n0, n1);
		}
	}
	(void) fputs("\n", output);

	(void) fputs("Residue Summary\n", output);
	for (p = rList_.first(); p != 0; rList_.next(p))
		rList_(p)->printSummary(output);
}

//
// Print pdb HELIX records
// id is a zero-based counter of the number of HELIX records printed
//
int
Model::printHelix(FILE *output, int id) const
{
	PDB pdb(PDB::HELIX);

	PDB::Helix &helix = pdb.helix;
	helix.comment[0] = '\0';
	for (Pix p = helixList_.first(); p != 0; helixList_.next(p)) {
		id++;
		Helix *h = helixList_(p);
		helix.serialNum = id;
		(void) sprintf(helix.id, "%d", id);
		helix.residues[0] = residue(h->from())->residue();
		helix.residues[1] = residue(h->to())->residue();
		helix.type = h->type();
		(void) fprintf(output, "%-71.71s%5d\n", pdb.chars(),
				h->to() - h->from() + 1);
	}
	return id;
}

//
// Print pdb SHEET records
// sid is a zero-based counter of the number of HELIX records printed
//
char
Model::printSheet(FILE *output, int sid) const
{
	//
	// Printing the sheet records is a bit tricky because
	// we need to derive the strands from the ladders.
	// We do so by first creating an ordered array of the
	// ladders.  We then define the first strand (method
	// depends on whether the sheet is cyclic or not), and
	// then iterate through.  We also have to do some post-
	// processing if the sheet is cyclic
	//
	PDB pdb(PDB::SHEET);
	PDB::Sheet &sheet = pdb.sheet;
	for (Pix p = sheetList_.first(); p != 0; sheetList_.next(p)) {
		Sheet *s = sheetList_(p);
		Ladder *fl = s->firstLadder();
		Ladder *pl = NULL;
		Ladder *spl = NULL;
		int ladderCount = s->ladderList().count();
		Ladder **lList = new Ladder *[ladderCount];
		int i = 0;
		for (Ladder *l = fl; l != NULL && !(l == fl && pl != NULL);
		spl = l, l = l->otherNeighbor(pl), pl = spl)
			lList[i++] = l;
		if (i != ladderCount) {
			(void) fprintf(stderr,
				"Inconsistent ladder count for sheet %c "
				"(%d should be %d)\n",
				s->name(), i, ladderCount);
			ladderCount = i;
		}
		int cyclic = fl->neighborCount() > 1;
		int overlap[2];

		PDB firstStrand(PDB::SHEET);
		PDB::Sheet &firstSheet = firstStrand.sheet;
		char sheetId[4];
		int idLen = 0;
		int idCount = sid++;
		sheetId[idLen++] = 'A' + (idCount % 26);
		idCount = idCount / 26;
		while (idLen < 3 && idCount > 0) {
			sheetId[idLen++] = 'A' + (idCount % 26) - 1;
			idCount = idCount / 26;
		}
		char *idPtr = firstSheet.id;
		while (--idLen >= 0)
			*idPtr++ = sheetId[idLen];
		*idPtr++ = '\0';
		firstSheet.count = ladderCount;
		firstSheet.strandNum = 1;
		firstSheet.sense = 0;
		if (cyclic) {
			pl = lList[ladderCount - 1];
			(void) fl->overlaps(pl, overlap);
			int start = min(fl->start(overlap[0]),
					pl->start(overlap[1]));
			int end = max(fl->end(overlap[0]),
					pl->end(overlap[1]));
			firstSheet.residues[0] = residue(start)->residue();
			firstSheet.residues[1] = residue(end)->residue();
			(void) fprintf(output, "%s\n", firstStrand.chars());
			registerLadder(pl, &firstSheet, overlap[1]);
		}
		else {
			firstSheet.count++;
			if (ladderCount == 1)
				// If there is only one ladder, then we print
				// strand 0 first, then strand 1
				overlap[0] = 0;
			else {
				// If there are more than one ladder, then
				// we find the overlapping strand with the
				// next ladder and use the other one
				(void) fl->overlaps(lList[1], overlap);
				overlap[0] = 1 - overlap[0];
			}
			firstSheet.residues[0] =
				residue(fl->start(overlap[0]))->residue();
			firstSheet.residues[1] =
				residue(fl->end(overlap[0]))->residue();
			(void) fprintf(output, "%s\n", firstStrand.chars());
		}

		sheet = firstSheet;
		for (i = 1; i < ladderCount; i++) {
			Ladder *l = lList[i];
			pl = lList[i - 1];
			(void) l->overlaps(pl, overlap);
			sheet.strandNum++;
			int start = min(l->start(overlap[0]),
					pl->start(overlap[1]));
			int end = max(l->end(overlap[0]),
					pl->end(overlap[1]));
			sheet.residues[0] = residue(start)->residue();
			sheet.residues[1] = residue(end)->residue();
			registerLadder(pl, &sheet, 1 - overlap[1]);
			(void) fprintf(output, "%s\n", pdb.chars());
		}

		if (cyclic)
			(void) fprintf(output, "%s\n", firstStrand.chars());
		else {
			sheet.strandNum++;
			int n = 1 - overlap[0];
			Ladder *&last = lList[ladderCount - 1];
			sheet.residues[0] = residue(last->start(n))->residue();
			sheet.residues[1] = residue(last->end(n))->residue();
			registerLadder(last, &sheet, overlap[0]);
			(void) fprintf(output, "%s\n", pdb.chars());
		}
		delete [] lList;
	}
	return sid;
}

//
// Set minimum number of residues in a helix
//
void
Model::setMinHelixLength(int n)
{
	minHelixLength = n;
}

//
// Set minimum number of residues in a strand
//
void
Model::setMinStrandLength(int n)
{
	minStrandLength = n;
}

//
// Set whether we should try to merge ladders in beta-bulges
//
void
Model::ignoreBulges(void)
{
	checkBulges = 0;
}

//
// Add the imide hydrogens to all residue
//
void
Model::addImideHydrogens(void)
{
	Pix p = rList_.first();
	Residue *prev = rList_(p);
	for (rList_.next(p); p != 0; rList_.next(p)) {
		Residue *r = rList_(p);
		(void) r->addImideHydrogen(prev);
		if (r->flag(R_TER))
			prev = NULL;
		else
			prev = r;
	}
}

//
// Find hydrogen bonds
//
void
Model::findHBonds(void)
{
	int max = rArray_->count();
	for (int i = 0; i < max; i++)
		for (int j = i + 2; j < max; j++) {
			(*hBond_)(i, j) = residue(i)->hBondedTo(residue(j));
			(*hBond_)(j, i) = residue(j)->hBondedTo(residue(i));
		}
}

//
// Find the n-turns (n = 3,4)
//
void
Model::findTurns(int n)
{
	int donor = n == 3 ? R_3DONOR : R_4DONOR;
	int acceptor = n == 3 ? R_3ACCEPTOR : R_4ACCEPTOR;
	int gap = n == 3 ? R_3GAP : R_4GAP;
	int max = rArray_->count() - n;
	for (int i = 0; i < max; i++)
		if (hBonded(i, i + n)) {
			residue(i)->setFlag(acceptor);
			for (int j = 1; j < n; j++)
				residue(i + j)->setFlag(gap);
			residue(i + n)->setFlag(donor);
		}
}

//
// Mark helices based on n-turn information
//
void
Model::markHelices(int n)
{
	int donor = n == 3 ? R_3DONOR : R_4DONOR;
	int acceptor = n == 3 ? R_3ACCEPTOR : R_4ACCEPTOR;
	int gap = n == 3 ? R_3GAP : R_4GAP;
	int helix = n == 3 ? R_3HELIX : R_4HELIX;
	int max = rArray_->count() - n;
	for (int i = 1; i < max; i++)
		if (residue(i - 1)->flag(acceptor)
		&&  residue(i)->flag(acceptor))
			for (int j = 0; j < n; j++)
				residue(i + j)->setFlag(helix);
}

//
// Construct helices based on marker information
//
void
Model::findHelices(void)
{
	int max = rArray_->count();
	int first = -1;
	for (int i = 0; i < max; i++)
		if (residue(i)->flag(R_3HELIX | R_4HELIX)) {
			if (first < 0)
				first = i;
		}
		else if (first >= 0) {
			if (i - first >= minHelixLength) {
				Helix *h = new Helix(first, i - 1);
				helixList_.append(h);
				h->setType(helixClass(h));
			}
			first = -1;
		}
}

//
// Find bridges
//
void
Model::findBridges(void)
{
	int max = rArray_->count();

	// First we construct a matrix and mark the bridges
	SquareArray<char> bridge(max);
	bridge.zero();
	int i;
	for (i = 1; i < max; i++) {
		for (int j = i + 1; j < max; j++) {
			if ((hBonded(i - 1, j) && hBonded(j, i + 1))
			||  (hBonded(j - 1, i) && hBonded(i, j + 1))) {
				bridge(i, j) = 'P';
				residue(i)->setFlag(R_PBRIDGE);
				residue(j)->setFlag(R_PBRIDGE);
			}
			else if ((hBonded(i, j) && hBonded(j, i))
			|| (hBonded(i - 1, j + 1) && hBonded(j - 1, i + 1))) {
				bridge(i, j) = 'A';
				residue(i)->setFlag(R_ABRIDGE);
				residue(j)->setFlag(R_ABRIDGE);
			}
		}
	}

	// Now we loop through and find the ladders
	int k;
	for (i = 0; i < max; i++) {
		for (int j = i + 1; j < max; j++) {
			switch (bridge(i, j)) {
			  case 'P':
				for (k = 0; bridge(i + k, j + k) == 'P'; k++) 
					bridge(i + k, j + k) = 'p';
				k--;
				ladderList_.append(new Ladder(B_PARA,
								i, i + k,
								j, j + k));
				break;
			  case 'A':
				for (k = 0; bridge(i + k, j - k) == 'A'; k++) 
					bridge(i + k, j - k) = 'a';
				k--;
				ladderList_.append(new Ladder(B_ANTI,
								i, i + k,
								j - k , j));
				break;
			}
		}
	}

	// Now we merge ladders of beta-bulges
	if (checkBulges)
		while (findBetaBulge())
			continue;

	// Finally we get rid of any ladder that is too short
	// (on either strand)
	int pruned;
	do {
		pruned = 0;
		for (Pix p = ladderList_.first(); p != 0; ladderList_.next(p)) {
			Ladder *l = ladderList_(p);
			if (l->end(0) - l->start(0) + 1 < minStrandLength
			||  l->end(1) - l->start(1) + 1 < minStrandLength) {
				ladderList_.remove(l);
				pruned = 1;
				break;
			}
		}
	} while (pruned);
}

//
// Find beta-bulges and merge the ladders
//
int
Model::findBetaBulge(void)
{
	for (Pix p1 = ladderList_.first(); p1 != 0; ladderList_.next(p1)) {
		Ladder *l1 = ladderList_(p1);
		if (l1->isBulge())
			continue;
		Pix p2 = p1;
		for (ladderList_.next(p2); p2 != 0; ladderList_.next(p2)) {
			Ladder *l2 = ladderList_(p2);
			if (l2->isBulge())
				continue;
			Ladder *l = Ladder::mergeBulge(l1, l2);
			if (l != NULL) {
				ladderList_.remove(l1);
				ladderList_.remove(l2);
				ladderList_.append(l);
				return 1;
			}
		}
	}
	return 0;
}

//
// Find beta-sheet based on ladder information
//
void
Model::findSheets(void)
{
	char sName = 'A';
	for (Pix p = ladderList_.first(); p != 0; ladderList_.next(p)) {
		Ladder *l = ladderList_(p);
		if (l->sheet() != NULL)
			continue;
		Sheet *s = new Sheet(sName);
		sheetList_.append(s);
		if (sName == 'Z')
			sName = 'A';
		else
			sName++;
		markLadder(l, s);
	}
}

//
// Mark this ladder (and all overlapped ladders, including closure)
// as part of given sheet
//
void
Model::markLadder(Ladder *ladder, Sheet *sheet)
{
	sheet->addLadder(ladder);
	ladder->setSheet(sheet);
	for (Pix p = ladderList_.first(); p != 0; ladderList_.next(p)) {
		Ladder *l = ladderList_(p);
		if (l->sheet() != NULL)
			continue;
		int overlap[2];
		if (!l->overlaps(ladder, overlap))
			continue;
		if (l->neighbor(overlap[0]) != NULL) {
			reportOverlap(l, overlap[0],
					ladder, l->neighbor(overlap[0]));
			continue;
		}
		if (ladder->neighbor(overlap[1]) != NULL) {
			reportOverlap(ladder, overlap[1],
					l, ladder->neighbor(overlap[1]));
			continue;
		}
		l->setNeighbor(overlap[0], ladder);
		ladder->setNeighbor(overlap[1], l);
		markLadder(l, sheet);
	}
}

//
// Report a ladder as overlapping on the same side with two other ladders
//
void
Model::reportOverlap(const Ladder *l, int side,
			const Ladder *o1, const Ladder *o2) const
{
	const PDB::Residue &first = residue(l->start(side))->residue();
	const PDB::Residue &last = residue(l->end(side))->residue();
	const PDB::Residue &ofirst = residue(l->start(1 - side))->residue();
	const PDB::Residue &olast = residue(l->end(1 - side))->residue();
	(void) fprintf(stderr,
		"Strand %d%c[%c]-%d%c[%c] (%d%c[%c]-%d%c[%c]) "
		"is paired with multiple ladders\n",
		first.seqNum, first.chainId, first.insertCode,
		last.seqNum, last.chainId, last.insertCode,
		ofirst.seqNum, ofirst.chainId, ofirst.insertCode,
		olast.seqNum, olast.chainId, olast.insertCode);

	const PDB::Residue &s10 = residue(o1->start(0))->residue();
	const PDB::Residue &e10 = residue(o1->end(0))->residue();
	const PDB::Residue &s11 = residue(o1->start(1))->residue();
	const PDB::Residue &e11 = residue(o1->end(1))->residue();
	(void) fprintf(stderr,
		"\t1 - Ladder %d%c[%c]-%d%c[%c], %d%c[%c]-%d%c[%c]\n",
		s10.seqNum, s10.chainId, s10.insertCode,
		e10.seqNum, e10.chainId, e10.insertCode,
		s11.seqNum, s11.chainId, s11.insertCode,
		e11.seqNum, e11.chainId, e11.insertCode);

	const PDB::Residue &s20 = residue(o2->start(0))->residue();
	const PDB::Residue &e20 = residue(o2->end(0))->residue();
	const PDB::Residue &s21 = residue(o2->start(1))->residue();
	const PDB::Residue &e21 = residue(o2->end(1))->residue();
	(void) fprintf(stderr,
		"\t2 - Ladder %d%c[%c]-%d%c[%c], %d%c[%c]-%d%c[%c]\n",
		s20.seqNum, s20.chainId, s20.insertCode,
		e20.seqNum, e20.chainId, e20.insertCode,
		s21.seqNum, s21.chainId, s21.insertCode,
		e21.seqNum, e21.chainId, e21.insertCode);
}

//
// Generate registration information for given ladder
//
void
Model::registerLadder(const Ladder *l, PDB::Sheet *sheet, int prev) const
{
	int cur = 1 - prev;
	if (l->type() == B_PARA) {
		//
		// We know that hBond(l->start(prev), l->start(cur))
		//
		sheet->sense = 1;
		Residue *r = residue(l->start(prev));
		if (r->hBondedTo(residue(l->start(cur) + 1))) {
			(void) strcpy(sheet->atoms[1].name, " O");
			sheet->atoms[1].residue = r->residue();
			r = residue(l->start(cur) + 1);
			(void) strcpy(sheet->atoms[0].name, " N");
			sheet->atoms[0].residue = r->residue();
		}
		else {
			r = residue(l->start(prev) + 1);
			(void) strcpy(sheet->atoms[1].name, " O");
			sheet->atoms[1].residue = r->residue();
			r = residue(l->start(cur));
			(void) strcpy(sheet->atoms[0].name, " N");
			sheet->atoms[0].residue = r->residue();
		}
	}
	else {
		//
		// We know that hBond(l->start(prev), l->end(cur))
		//
		sheet->sense = -1;
		Residue *r = residue(l->start(prev));
		if (r->hBondedTo(residue(l->end(cur)))) {
			(void) strcpy(sheet->atoms[1].name, " O");
			sheet->atoms[1].residue = r->residue();
			r = residue(l->end(cur));
			(void) strcpy(sheet->atoms[0].name, " N");
			sheet->atoms[0].residue = r->residue();
		}
		else {
			r = residue(l->start(prev) + 1);
			(void) strcpy(sheet->atoms[1].name, " O");
			sheet->atoms[1].residue = r->residue();
			r = residue(l->end(cur) - 1);
			(void) strcpy(sheet->atoms[0].name, " N");
			sheet->atoms[0].residue = r->residue();
		}
	}
}

//
// Determine the PDB class of the helix
//
int
Model::helixClass(const Helix *h) const
{
	Atom *ca[4];
	int from = h->from();
	Residue *r = residue(from);
	for (int i = 0; i < 4; i++) {
		ca[i] = residue(from + i)->atom(" CA");
		if (ca[i] == NULL)
			return 0;
	}
	float angle = dihedral(ca[0]->coord(), ca[1]->coord(),
				ca[2]->coord(), ca[3]->coord());
	if (angle > 0) {
		if (r->flag(R_4HELIX))
			return 1;
		else if (r->flag(R_3HELIX))
			return 5;
	}
	else {
		if (r->flag(R_4HELIX))
			return 6;
	}
	return 0;
}
