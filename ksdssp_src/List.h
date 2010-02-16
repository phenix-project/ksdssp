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

#ifndef LIST_INCLUDE
#define	LIST_INCLUDE

#include <stddef.h>

typedef void	*Pix;

template <class T>
class ListItem	{
	ListItem<T>	*_next;
	T		*_value;
public:
			ListItem(T *value)
				{
					_next = NULL;
					_value = value;
				}
	void		append(ListItem<T> *element)
				{
					ListItem<T> *n;
					for (n = this; n->_next != NULL;
					n = n->_next)
						continue;
					n->_next = element;
				}
	ListItem<T>	*next(void) const
				{
					return _next;
				}
	void		setNext(ListItem<T> *n)
				{
					_next = n;
				}
	T		*value(void) const
				{
					return _value;
				}
};

template <class T>
class List	{
	ListItem<T>	*_head;
public:
			List(void);
			List(List *orig);
			~List(void);
	int		add(T *element);
	int		append(T *element);
	int		remove(T *element);
	void		clear(void);
	T		*head(void);
	T		*pop(void);
	T		*find(const void *data, int (*match)(const T *t,
						const void *mdata)) const;

	int		count(void) const;
	Pix		first(void) const;
	void		next(Pix &p) const;
	T		*operator()(Pix p) const;
};

//
// Implementation of List
//

template <class T>
List<T>::List(void)
{
	_head = NULL;
}

template <class T>
List<T>::List(List *orig)
{
	// This constructor is mainly used in functions that return
	// a list.  The function can create a local list, then
	//	return List(&local_list);
	// which should prevent a deep copy (that we don't currently provide)
	_head = orig->_head;
	orig->_head = NULL;
}

template <class T>
List<T>::~List(void)
{
	clear();
}

template <class T>
int
List<T>::add(T *element)
{
	ListItem<T> *item = new ListItem<T>(element);
	item->append(_head);
	_head = item;
	return 0;
}

template <class T>
int
List<T>::append(T *element)
{
	ListItem<T> *item = new ListItem<T>(element);
	if (_head == NULL)
		_head = item;
	else
		_head->append(item);
	return 0;
}

template <class T>
int
List<T>::remove(T *element)
{
	ListItem<T> *p = NULL;
	ListItem<T> *t;
	for (t = _head; t != NULL; p = t, t = t->next())
		if (t->value() == element)
			break;
	if (t == NULL)
		return -1;
	if (p == NULL)
		_head = t->next();
	else
		p->setNext(t->next());
	delete t;
	return 0;
}

template <class T>
void
List<T>::clear(void)
{

	while (_head != NULL) {
		ListItem<T> *n = _head->next();
		delete _head;
		_head = n;
	}
}

template <class T>
T *
List<T>::head(void)
{
	if (_head == NULL)
		return NULL;
	return _head->value();
}

template <class T>
T *
List<T>::pop(void)
{
	if (_head == NULL)
		return NULL;
	T *v = _head->value();
	ListItem<T> *n = _head->next();
	delete _head;
	_head = n;
	return v;
}

template <class T>
T *
List<T>::find(const void *data, int (*match)(const T *t, const void *mdata))
		const
{
	for (Pix p = first(); p != 0; next(p)) {
		T *element = (*this)(p);
		if ((*match)(element, data))
			return element;
	}
	return NULL;
}

template <class T>
int
List<T>::count(void) const
{
	int count = 0;
	for (ListItem<T> *item = _head; item != NULL; item = item->next())
		count++;
	return count;
}

template <class T>
Pix	
List<T>::first(void) const
{
	return _head;
}

template <class T>
void
List<T>::next(Pix &p) const
{
	p = ((ListItem<T> *) p)->next();
}

template <class T>
T *
List<T>::operator()(Pix p) const
{
	return ((ListItem<T> *) p)->value();
}

#endif
