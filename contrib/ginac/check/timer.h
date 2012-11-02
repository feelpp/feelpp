/** @file timer.h
 *
 *  A simple stop watch class. */

/*
 *  GiNaC Copyright (C) 1999-2011 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef TIMER_H
#define TIMER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_RUSAGE
#include <sys/resource.h>
#else
#include <ctime>
#endif

class timer {
public:
	timer();
	void start();
	void stop();
	void reset();
	double read();
	bool running();
private:
	bool on;
#ifdef HAVE_RUSAGE
	struct rusage used1, used2;
#else
	std::clock_t used1, used2;
#endif
};

#endif // ndef TIMER_H
