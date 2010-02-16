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

#ifndef listarray_h
#define listarray_h

#include "List.h"

template <class T>
class ListArray {
	int		size_;
	T		**array_;
public:
			ListArray(const List<T> &list);
			~ListArray(void);
	T		*operator()(int n) const;
	int		count(void) const;
};

template <class T>
ListArray<T>::ListArray(const List<T> &list)
{
	size_ = list.count();
	array_ = new T *[size_];
	int i = 0;
	for (Pix p = list.first(); p != 0; list.next(p))
		array_[i++] = list(p);
}

template <class T>
ListArray<T>::~ListArray(void)
{
	delete[] array_;
}

template <class T>
T *
ListArray<T>::operator()(int n) const
{
	if (n < 0 || n >= size_)
		return NULL;
	return array_[n];
}

template <class T>
int
ListArray<T>::count(void) const
{
	return size_;
}

#endif
