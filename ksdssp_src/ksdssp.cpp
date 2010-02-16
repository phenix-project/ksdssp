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

#ifdef NEED_BOOL
#define	bool	int
#define	false	0
#define	true	1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <pdb++.h>
#include "ksdssp.h"
#include "Model.h"
#include "XGetopt.h"

#ifndef DONT_INSTANIATE
template class List<Model>;
#endif

#if defined(NeXT) || defined(mips)
extern "C" char *strerror(int);
#endif

int verbose = 0;

//
// This is an implementation of
//
//	Dictionary of Protein Secondary Structure:
//	Pattern Recognition of Hydrogen-Bonded and
//	Geometrical Features
//	Wolfgang Kabsch and Christian Sander
//	Biopolymers, Vol. 22, 2577-2637 (1983)
//
int
main(int argc, char **argv)
{
	// Parse command line options
	int o;
	char *summaryFile = NULL;
	while ((o = Xgetopt(argc, argv, "c:h:s:vBS:")) != EOF)
		switch (o) {
		  case 'c':
			Residue::setHBondCutoff(atof(optarg));
			break;
		  case 'h':
			Model::setMinHelixLength(atoi(optarg));
			break;
		  case 's':
			Model::setMinStrandLength(atoi(optarg));
			break;
		  case 'v':
			verbose++;
			break;
		  case 'B':
			Model::ignoreBulges();
			break;
		  case 'S':
			summaryFile = optarg;
			break;
		}

	// Check input PDB file
	FILE *input = NULL;
	FILE *output = NULL;
	const char *inputFile = NULL;
	const char *outputFile = NULL;
	switch (argc - optind) {
	  case 0:
		input = stdin;
		inputFile = "standard input";
		output = stdout;
		outputFile = "standard output";
		break;
	  case 1:
		inputFile = argv[optind];
		if (strcmp(inputFile, "-") == 0) {
			input = stdin;
			inputFile = "standard input";
		}
		else {
			input = fopen(inputFile, "r");
			if (input == NULL) {
				(void) fprintf(stderr, "%s: %s: %s\n",
					argv[0], inputFile, strerror(errno));
				return 1;
			}
		}
		output = stdout;
		outputFile = "standard output";
		break;
	  case 2:
		inputFile = argv[optind];
		if (strcmp(inputFile, "-") == 0) {
			input = stdin;
			inputFile = "standard input";
		}
		else {
			input = fopen(inputFile, "r");
			if (input == NULL) {
				(void) fprintf(stderr, "%s: %s: %s\n",
					argv[0], argv[optind],
					strerror(errno));
				return 1;
			}
		}
		outputFile = argv[optind + 1];
		if (strcmp(outputFile, "-") == 0) {
			output = stdout;
			outputFile = "standard output";
		}
		else {
			output = fopen(outputFile, "w");
			if (output == NULL) {
				(void) fprintf(stderr, "%s: %s: %s\n",
					argv[0], outputFile, strerror(errno));
				return 1;
			}
		}
		break;
	  default:
		(void) fprintf(stderr, "Usage: %s [-c config] [pdb_file]\n",
			argv[0]);
		return 1;
	}

	// Construct molecule from PDB file
	List<Model> modelList;
	for (;;) {
		Model *m = new Model(input);
		if (!m->okay()) {
			(void) fprintf(stderr, "%s: %s: %s\n",
				argv[0], inputFile, m->error());
			return 1;
		}
		int anyMore = m->anyMore();
		if (m->anyAtoms())
			modelList.append(m);
		else
			delete m;
		if (!anyMore)
			break;
	}
	if (modelList.count() <= 0) {
		(void) fprintf(stderr, "%s: %s: no atoms read\n",
			argv[0], inputFile);
		return 1;
	}

	// Compute secondary structure and print helix and sheet records
	Pix p;
	for (p = modelList.first(); p != 0; modelList.next(p))
		modelList(p)->defineSecondaryStructure();
	Pix hp = modelList.first();
	Pix sp = hp;
	int fileCount = 0;
	while (hp != 0) {
		if (fileCount++ > 0)
			(void) fprintf(output, "%s\n", PDB(PDB::END).chars());
		int helixId = 0;
		int sheetId = 0;
		Model *m = modelList(hp);
		if (m->modelNumber() != -1)
			fprintf(output, "%s\n", m->fileRecord().chars());
		int modelNumber = m->modelNumber();
		for (; hp != 0 && modelList(hp)->modelNumber() == modelNumber;
		modelList.next(hp))
			helixId = modelList(hp)->printHelix(output, helixId);
		for (; sp != hp; modelList.next(sp))
			sheetId = modelList(sp)->printSheet(output, sheetId);
	}
	if (fileCount > 1)
		(void) fprintf(output, "%s\n", PDB(PDB::END).chars());

	// Print chain summaries
	FILE *summary = NULL;
	if (summaryFile != NULL
	&& (summary = fopen(summaryFile, "w")) == NULL) {
		(void) fprintf(stderr, "%s: %s: %s\n",
			argv[0], summaryFile, strerror(errno));
		return 1;
	}
	if (summary != NULL) {
		for (p = modelList.first(); p != 0; modelList.next(p))
			modelList(p)->printSummary(summary);
		(void) fclose(summary);
	}

	return 0;
}
