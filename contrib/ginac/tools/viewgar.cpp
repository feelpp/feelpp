/** @file viewgar.cpp
 *
 *  GiNaC archive file viewer. */

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

#include "ginac.h"
using namespace GiNaC;

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " [-d] file..." << endl;
		exit(1);
	}
	--argc; ++argv;

	bool dump_mode = false;
	try {
		lst l;
		while (argc) {
			if (strcmp(*argv, "-d") == 0) {
				dump_mode = true;
				--argc; ++argv;
			}
			std::ifstream f(*argv, std::ios_base::binary);
			archive ar;
			f >> ar;
			if (dump_mode) {
				ar.printraw(std::cout);
				std::cout << std::endl;
			} else {
				for (unsigned int i=0; i<ar.num_expressions(); ++i) {
					std::string name;
					ex e = ar.unarchive_ex(l, name, i);
					std::cout << name << " = " << e << std::endl;
				}
			}
			--argc; ++argv;
		}
	} catch (std::exception &e) {
		std::cerr << *argv << ": " << e.what() << std::endl;
	}
}
