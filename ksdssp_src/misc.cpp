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
#include "ksdssp.h"
#include "misc.h"

//
// Compute distance squared between two points
//
float
distSquared(const float v1[3], const float v2[3])
{
	float l = 0;
	for (int i = 0; i < 3; i++) {
		float d = v2[i] - v1[i];
		l += d * d;
	}
	return l;
}

//
// Compute distance between two points
//
float
distance(const float v1[3], const float v2[3])
{
	float l = 0;
	for (int i = 0; i < 3; i++) {
		float d = v2[i] - v1[i];
		l += d * d;
	}
	return sqrtf(l);
}

//
// Compute distance between two points
//
float
magnitude(const float v[3])
{
	float l = 0;
	for (int i = 0; i < 3; i++)
		l += v[i] * v[i];
	return sqrtf(l);
}

//
// Normalize given vector
//
void
normalize(float r[3])
{
	float l = magnitude(r);
	for (int i = 0; i < 3; i++)
		r[i] /= l;
}

//
// Bisect two unit vectors
//
void
bisect(float r[3], const float v1[3], const float v2[3])
{
	for (int i = 0; i < 3; i++)
		r[i] = (v1[i] + v2[i]) / 2;
	normalize(r);
}

//
// Compute angle between two vectors
//
float
angle(const float v1[3], const float v2[3])
{
	return acos(dotProduct(v1, v2) / magnitude(v1) / magnitude(v2));
}

//
// Compute angle defined by three points
//
float
angle(const float v1[3], const float v2[3], const float v3[3])
{
	float d1[3], d2[3];
	for (int i = 0; i < 3; i++) {
		d1[i] = v1[i] - v2[i];
		d2[i] = v3[i] - v2[i];
	}
	return angle(d1, d2);
}

//
// Compute the cross product of two vectors
//
void
crossProduct(float answer[3], const float v1[3], const float v2[3])
{
	answer[0] = v1[1] * v2[2] - v2[1] * v1[2];
	answer[1] = v1[2] * v2[0] - v2[2] * v1[0];
	answer[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

//
// Compute the dot product of two vectors
//
float
dotProduct(const float v1[3], const float v2[3])
{
	float d = 0;
	for (int i = 0; i < 3; i++)
		d += v1[i] * v2[i];
	return d;
}

//
// Compute the dihedral formed by four vectors
//
float
dihedral(const float v1[3], const float v2[3],
		const float v3[3], const float v4[3])
{
	float d12[3], d32[3];
	float d23[3], d43[3];
	for (int i = 0; i < 3; i++) {
		d12[i] = v1[i] - v2[i];
		d32[i] = v3[i] - v2[i];
		d23[i] = v2[i] - v3[i];
		d43[i] = v4[i] - v3[i];
	}

	float d1[3], d2[3];
	crossProduct(d1, d12, d32);
	crossProduct(d2, d23, d43);
	return angle(d1, d2);
}
