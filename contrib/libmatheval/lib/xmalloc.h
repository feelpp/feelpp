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

#ifndef XMALLOC_H
#define XMALLOC_H 1

/* Macro definitions to simplify corresponding function calls.  */
#define XMALLOC(type, num) ((type *) xmalloc ((num) * sizeof(type)))
#define XREALLOC(type, ptr, num) ((type *) xrealloc ((ptr), (num) * sizeof(type)))
#define XCALLOC(type, num) ((type *) xcalloc ((num), sizeof(type)))
#define XFREE(stale) free (stale);

/* Replacement for malloc() function with error checking.  */
void           *xmalloc(size_t size);

/* Same as above from realloc().  */
void           *xrealloc(void *ptr, size_t size);

/* Same as above for calloc().  */
void           *xcalloc(size_t num, size_t size);

#endif
