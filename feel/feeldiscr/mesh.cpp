/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2012-12-02

  Copyright (C) 2012 Université de Strasbourg

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
/**
   \file mesh.cpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2012-12-02
 */

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>


#include <map>

#include <iostream>
#include <string>
#include <complex>

#include <feel/feeldiscr/mesh.hpp>


namespace Feel
{
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;
}

// We need to tell fusion about our Model struct
// to make it a first-class fusion citizen. This has to
// be in global scope.

BOOST_FUSION_ADAPT_STRUCT(
	Feel::MeshMarkerName,
	(std::string, name)
	(std::vector<int>, ids)
	)
//]

namespace Feel
{
///////////////////////////////////////////////////////////////////////////////
//  Our Model parser
///////////////////////////////////////////////////////////////////////////////
//[tutorial_Model_parser
template <typename Iterator>
struct MarkerParser : qi::grammar<Iterator, std::vector<MeshMarkerName>(), ascii::space_type>
{
        MarkerParser() : MarkerParser::base_type(start,"start")
		{
			using qi::int_;
			using qi::lit;
			using ascii::char_;

			vints %= int_ % ',';
			marker %=  +(char_-'=')  >> '=' >> '{' >>  vints   >> '}';
			start %= '{' >> marker % ';' >>  '}' ;

		}
	qi::rule<Iterator, std::vector<int>(), ascii::space_type> vints;
	qi::rule<Iterator, MeshMarkerName(), ascii::space_type> marker;
	qi::rule<Iterator, std::vector<MeshMarkerName>(), ascii::space_type> start;
};
//]


std::vector<MeshMarkerName> markerMap( int Dim )
{
    using boost::spirit::ascii::space;
    typedef std::string::const_iterator iterator_type;
    typedef Feel::MarkerParser<iterator_type> MarkerParser;
    std::string section = (boost::format("mesh%1%d.markers")% Dim).str();

    MarkerParser g;
    std::vector<Feel::MeshMarkerName> emp;
    if ( Environment::vm().count(section) == 0 )
        return emp;
    std::string str = Environment::vm( _name=section ).as<std::string>();
    LOG(INFO) << "markers=" << str << "\n";

    std::string::const_iterator iter = str.begin();
    std::string::const_iterator end = str.end();
    bool r = phrase_parse(iter, end, g, space, emp);
    CHECK( r && ( iter == end) ) << "parsing mesh marker names failed\n";

    LOG(INFO) << emp.size() << "\n";
    for(int i = 0; i < emp.size(); ++i )
    {
        LOG(INFO) << "id {" << emp[i].ids[0] << "," <<  emp[i].ids[1] << "} marker: "<< emp[i].name << "\n";
    }
    google::FlushLogFiles(google::GLOG_INFO);
    return emp;
}

po::options_description mesh_options( int Dim, std::string const& prefix )
{
    po::options_description _options( (boost::format("Mesh %1%D ")%Dim).str() + prefix + " options" );
    _options.add_options()
        ( prefixvm( prefix, (boost::format("mesh%1%d.hsize") % Dim).str() ).c_str(), Feel::po::value<double>()->default_value(0.1), "mesh characteristic size" )
        ( prefixvm( prefix, (boost::format("mesh%1%d.markers") % Dim).str() ).c_str(), Feel::po::value<std::string>(), "mesh markers map" )
        ( prefixvm( prefix, (boost::format("mesh%1%d.localisation.use-extrapolation") % Dim).str() ).c_str(), Feel::po::value<bool>()->default_value(true),
          "use extrapolation if localisation fails" )
        ( prefixvm( prefix, (boost::format("mesh%1%d.localisation.nelt-in-leaf-kdtree") % Dim).str() ).c_str(), Feel::po::value<int>()->default_value(-1),
          "use extrapolation if localisation fails" )
        ;
    return _options;
}
} // Feel
