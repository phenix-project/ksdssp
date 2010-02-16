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

#ifndef squarearray_h
#define squarearray_h

#include <string.h>

template <class T>
class SquareArray
{
	int		dim_;
	T		*array_;
public:
			SquareArray(int size);
			~SquareArray(void);
	void		zero(void);
	int		dimension(void) const;
	T		&operator()(int row, int col);
};

template <class T>
SquareArray<T>::SquareArray(int size)
{
	dim_ = size;
	array_ = new T[dim_ * dim_];
}

template <class T>
SquareArray<T>::~SquareArray(void)
{
	delete[] array_;
}

template <class T>
int
SquareArray<T>::dimension(void) const
{
	return dim_;
}

template <class T>
void
SquareArray<T>::zero(void)
{
	memset(array_, 0, dim_ * dim_ * sizeof (T));
}

template <class T>
T &
SquareArray<T>::operator()(int row, int col)
{
	return array_[row * dim_ + col];
}

#endif
