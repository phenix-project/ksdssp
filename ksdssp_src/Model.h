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

#ifndef model_h
#define model_h

#include <stdio.h>
#include <string>
#include "Residue.h"
#include "List.h"
#include "ListArray.h"
#include "SquareArray.h"
#include "Structure.h"

class Model {
	int			anyMore_;
	std::string		error_;
	List<Residue>		rList_;
	ListArray<Residue>	*rArray_;
	SquareArray<char>	*hBond_;
	List<Helix>		helixList_;
	List<Ladder>		ladderList_;
	List<Sheet>		sheetList_;
	int			modelNumber_;
	PDB			fileRecord_;
public:
			Model(FILE *input);
			~Model(void);
	int		okay(void) const { return error_ == ""; }
	int		anyMore(void) const { return anyMore_; }
	int		anyAtoms(void) const { return rList_.count() > 0; }
	const char	*error(void) const { return error_.c_str(); }
	int		modelNumber(void) const { return modelNumber_; }
	const PDB	&fileRecord(void) const { return fileRecord_; }
	void		defineSecondaryStructure(void);
	void		printResidues(FILE *output) const;
	void		printSummary(FILE *output) const;
	int		printHelix(FILE *output, int id) const;
	char		printSheet(FILE *output, int id) const;
public:
	static void	setMinStrandLength(int n);
	static void	setMinHelixLength(int n);
	static void	ignoreBulges(void);
private:
	int		hBonded(int i, int j) { return (*hBond_)(i, j); }
	Residue		*residue(int n) const { return (*rArray_)(n); }
	void		addImideHydrogens(void);
	void		findHBonds(void);
	void		findTurns(int n);
	void		markHelices(int n);
	void		findHelices(void);
	void		findBridges(void);
	int		findBetaBulge(void);
	void		findSheets(void);
	void		markLadder(Ladder *ladder, Sheet *sheet);
	void		reportOverlap(const Ladder *l, int s, const Ladder *o1,
					const Ladder *o2) const;
	void		registerLadder(const Ladder *l, PDB::Sheet *sheet,
					int prev) const;
	int		helixClass(const Helix *h) const;
};

#endif
