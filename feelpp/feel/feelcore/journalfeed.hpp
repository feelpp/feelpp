// This file is part of the Feel library
//
// Author(s): Feel++ Contortium
//      Date: 2017-07-10
//
// @copyright (C) 2017 University of Strasbourg
// @copyright (C) 2012-2017 Feel++ Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef FEELPP_JOURNALFEED_HPP
#define FEELPP_JOURNALFEED_HPP 1

#include <map>
#include <initializer_list>
#include <boost/variant.hpp> // C++17 to std.
#include <boost/lexical_cast.hpp>
#include <feel/feelcore/journalwatcher.hpp>

// All Feel++ observers are defined here.
namespace Feel
{

// Feed the journal with map entries.
// Example:
//     To add two entries a, b in the ptree and Json.
//     JournalFeel{{"a",1},{"b,2"}};
class JournalFeed : virtual public JournalWatcher
{
    using super = JournalWatcher;
public:
    using value_type = boost::variant<std::string,int, double,
          std::vector<int>, 
          std::vector<double>,
          std::vector<std::string>>;
    using init_type = std::initializer_list< std::pair<std::string const, value_type> >;
    using map_type = std::map< std::string /*const*/, value_type >;

	//! Constructor.
    //! @{
    JournalFeed() : super("") {}

    //! Create a ptree from initializer list.
    explicit JournalFeed( init_type const& m )
        :
        JournalFeed( map_type(m) )
        {}

    //! Create a ptree from map.
    explicit JournalFeed( map_type const& m )
        :
        super("")
    {
        pt::ptree pt;
        for( const auto& kv: m )
        {
            // Visit map value (variant) and transform the value into string.
            const std::string s = boost::apply_visitor( JournalFeedVisitor(), kv.second );
			pt.put( kv.first, s );
        }
        this->setInformationObject( pt );
    }

    // @}

    void add( std::string const& key, value_type val )
        {
            const std::string s = boost::apply_visitor( JournalFeedVisitor(), val );
            pt::ptree pt;
            pt.put( key, s );
            this->addInformationObject( pt );
        }

    //! Return a pointer to this object.
    const JournalFeed* operator()()
    {
        return this;
    }

private:
	//! Visitor class which converts variant contained type to a string.
	struct JournalFeedVisitor : public boost::static_visitor<std::string>
	{
		// Handle containers into strings.
		template<typename T0, 
                 template<class T1, class = std::allocator<T1> >
                     class container_type = std::vector >
			const std::string operator()(const container_type<T0>& arg) const
			{
				if( not arg.empty() )
				{
					std::ostringstream oss;
					std::copy(arg.begin(), arg.end()-1,
							  std::ostream_iterator<T0>(oss, ","));
					oss << arg.back();
					return oss.str();
				}
				return "";
			}

		// Default
		template<typename T>
			const std::string operator()(T const& arg ) const
			{
				return boost::lexical_cast<std::string>(arg);
			}
	};

}; // JournalFeed class.

} // Feel namespace.

#endif // FEELPP_JOURNALFEED_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
