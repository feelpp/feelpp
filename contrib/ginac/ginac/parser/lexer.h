/** @file lexer.h
 *
 *  The lexer. */

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

#ifndef GINAC_LEXER_H
#define GINAC_LEXER_H

#include <iosfwd>
#include <string>
#include <cstddef>

namespace GiNaC {

class lexer
{
	std::istream* input;
	std::ostream* output;
	std::ostream* error;
	/// last character read from stream
	int c;
	/// identifier and number tokens are stored here
	std::string str;
	std::size_t line_num;
	std::size_t column;
	friend class parser;
public:

	lexer(std::istream* in = 0, std::ostream* out = 0, std::ostream* err = 0);
	~lexer();

	int gettok();
	void switch_input(std::istream* in);

	struct token_type
	{
		enum
		{
			eof		= -1,
			identifier	= -4,
			number		= -5,
			literal		= -6

		};
	};
	
	/// Symbolic name of the token (for error reporting)
	std::string tok2str(const int tok) const;
};

} // namespace GiNaC

#endif // ndef GINAC_LEXER_H
