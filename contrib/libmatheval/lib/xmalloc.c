/* 
 * Copyright (C) 1999, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2011
 * Free Software Foundation, Inc.
 * 
 * This file is part of GNU libmatheval
 * 
 * GNU libmatheval is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * GNU libmatheval is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU libmatheval.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "error.h"

#if !HAVE_BZERO && HAVE_MEMSET
#define bzero(buf, bytes) ((void) memset (buf, 0, bytes))
#endif

void           *
xmalloc(size_t num)
{
	/* Call malloc() and check return value. */
	void           *ptr_new = malloc(num);

	if (!ptr_new)
		error_fatal("unable to allocate memory");
	return ptr_new;
}

void           *
xrealloc(void *ptr, size_t num)
{
	void           *ptr_new;

	/* If memory not already allocated, fall back to xmalloc(). */
	if (!ptr)
		return xmalloc(num);

	/* Otherwise, call realloc() and check return value. */
	ptr_new = realloc(ptr, num);
	if (!ptr_new)
		error_fatal("unable to allocate memory");

	return ptr_new;
}

void           *
xcalloc(size_t num, size_t size)
{
	/* Allocate memory and fill with zeroes. */
	void           *ptr_new = xmalloc(num * size);

	bzero(ptr_new, num * size);
	return ptr_new;
}
