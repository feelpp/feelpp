/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 sept. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <feel/feelcore/removecomments.hpp>

#include <feel/feelcore/feel.hpp>

//  Copyright (c) 2001-2010 Hartmut Kaiser
//  Copyright (c) 2001-2007 Joel de Guzman
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/lex_lexertl.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>

#include <iostream>
#include <string>

namespace Feel
{

namespace spirit = boost::spirit;

///////////////////////////////////////////////////////////////////////////////
//  Token definition: We use the lexertl based lexer engine as the underlying
//                    lexer type.
///////////////////////////////////////////////////////////////////////////////
enum tokenids
{
    IDANY = spirit::lex::min_token_id + 10
};

template <typename Lexer>
struct strip_comments_tokens : spirit::lex::lexer<Lexer>
{
    strip_comments_tokens()
        : strip_comments_tokens::base_type( spirit::lex::match_flags::match_default )
    {
        // define tokens and associate them with the lexer
        cppcomment = "\"//\"[^\n]*"; // '//[^\n]*'
        ccomment = "\"/*\"";         // '/*'
        endcomment = "\"*/\"";       // '*/'

        // The following tokens are associated with the default lexer state
        // (the "INITIAL" state). Specifying 'INITIAL' as a lexer state is
        // strictly optional.
        this->self.add( cppcomment )   // no explicit token id is associated
            ( ccomment )( ".", IDANY ) // IDANY is the token id associated with this token
            // definition
            ;

        // The following tokens are associated with the lexer state "COMMENT".
        // We switch lexer states from inside the parsing process using the
        // in_state("COMMENT")[] parser component as shown below.
        this->self( "COMMENT" ).add( endcomment )( ".", IDANY );
    }

    spirit::lex::token_def<> cppcomment, ccomment, endcomment;
};

///////////////////////////////////////////////////////////////////////////////
//  Grammar definition
///////////////////////////////////////////////////////////////////////////////
template <typename Iterator>
struct strip_comments_grammar : spirit::qi::grammar<Iterator>
{
    template <typename TokenDef>
    strip_comments_grammar( TokenDef const& tok, std::ostringstream& ostr )
        : strip_comments_grammar::base_type( start )
    {
        using namespace spirit;
        // The in_state("COMMENT")[...] parser component switches the lexer
        // state to be 'COMMENT' during the matching of the embedded parser.
        start = *( tok.ccomment >> spirit::qi::in_state( "COMMENT" )
                                       [
                                           // the lexer is in the 'COMMENT' state during
                                           // matching of the following parser components
                                           *token( IDANY ) >> tok.endcomment] |
                   tok.cppcomment | qi::token( IDANY )[ostr << _1] );
    }

    spirit::qi::rule<Iterator> start;
};

std::string
removeComments( std::string str )
{
    using namespace spirit;
    // iterator type used to expose the underlying input stream
    typedef std::string::iterator base_iterator_type;

    // lexer type
    typedef lex::lexertl::lexer<lex::lexertl::token<base_iterator_type>>
        lexer_type;

    // iterator type exposed by the lexer
    typedef strip_comments_tokens<lexer_type>::iterator_type iterator_type;

    // now we use the types defined above to create the lexer and grammar
    // object instances needed to invoke the parsing process
    strip_comments_tokens<lexer_type> strip_comments;
    std::ostringstream ostr;
    strip_comments_grammar<iterator_type> g( strip_comments, ostr );

    // Parsing is done based on the token stream, not the character
    // stream read from the input.
    base_iterator_type first = str.begin();

    bool r = lex::tokenize_and_parse( first, str.end(), strip_comments, g );

    if ( r )
    {
        //std::cout << "removecomments parsing succeeded\n";
    }
    else
    {
        std::string rest( first, str.end() );
        LOG( INFO ) << "removecomments parsing failed\n";
        LOG( INFO ) << "stopped at: \"" << rest << "\"\n";
    }

    return ostr.str();
}

} // Feel
