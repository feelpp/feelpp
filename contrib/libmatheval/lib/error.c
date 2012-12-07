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

const char      lib_name[] = "libmatheval";	/* Library name to be
						 * printed to standard
						 * error stream in order
						 * to locate problem
						 * source.  */

/* Issue error message to standard error stream and, if status value not
 * less than zero, terminate calling program.  Argument mode is intended
 * to describe error severity. */
static void     error_issue(int status, const char *mode,
			    const char *message);

void
error_warning(const char *message)
{
	/* Issue warning. */
	error_issue(-1, "warning", message);
}

void
error_fatal(const char *message)
{
	/* Issue error. */
	error_issue(EXIT_FAILURE, "FATAL", message);
}

static void
error_issue(int status, const char *mode, const char *message)
{
	/* Print error message. */
	fprintf(stderr, "%s: %s: %s\n", lib_name, mode, message);

	/* Check status and eventually terminate calling program. */
	if (status >= 0)
		exit(status);
}
