/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 19 Feb 2016

 Copyright (C) 2016 Feel++ Consortium

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
#ifndef FEELPP_FEELIO_HPP
#define FEELPP_FEELIO_HPP 1

#include <iostream>
#include <fstream>
#include <feel/feelcore/environment.hpp>


namespace Feel {

/**
 * Output Stream that outputs only on master rank of a worldcomm
 *
 * The first application is to use the instantiation of MasterStream cout, cerr
 * or clog to default output on master rank process of a Feel++ application
 * @code
 * Environment env(...); // initialize Feel++ environment
 * // cout only on master rank process
 * cout << "Hello World from process " << Environment::rank()  << std::endl;
 * @encode
 */
class FEELPP_EXPORT MasterStream
{
public:
    /**
     * Construct a MasterStream from a std::ostream and a WorldComm. The
     * WorldComm master rank process is in charge of outputting the data from
     * the MasterStream
     */
    MasterStream(std::ostream& _out, std::shared_ptr<WorldComm> _wc = Environment::worldCommPtr() )
        :
        out(_out),
        wc(_wc)
        {}

    void attachWorldComm( std::shared_ptr<WorldComm>& _wc )
        {
            wc = _wc;
        }

    /**
     * this overload allows to output stream to the master rank process
     */
    template<typename T>
    const MasterStream& operator<<(const T& v) const
        {
            if ( wc->isMasterRank() )
                out << v;
            return *this;
        }
    /**
     * this overload handles std::endl and std::flush and allows to honor them
     */
    MasterStream const& operator<<(std::ostream& (*F)(std::ostream&)) const
        {
            if ( wc->isMasterRank() )
                F(out);
            return *this;
        }

    /**
     * provide an interface to \c str() for ostringstream
     *
     * @return the context of the string stream if it is an ostringstream, an
     * empty string otherwise
     */
    std::string str()
        {
            auto* o = dynamic_cast<std::ostringstream*>( &out );
            if ( o )
                return o->str();
            return std::string();
        }

    /**
     * @return the shared_ptr WorldComm
     */
    std::shared_ptr<WorldComm> worldCommPtr() const { return wc; }

    /**
     * variadic template function to provide open() from std::ofstream
     */
    template<typename... Args>
    void open(Args... args)
        {
            auto* o = dynamic_cast<std::ofstream*>( &out );
            if ( o && wc->isMasterRank() )
                o->open( args... );
        }

    /**
     * provides close() from std::ofstream
     */
    void close()
        {
            auto* o = dynamic_cast<std::ofstream*>( &out );
            if ( o && wc->isMasterRank() )
                o->close();
        }

protected:
    std::ostream& out;
    std::shared_ptr<WorldComm> wc;
};


extern MasterStream cout;
extern MasterStream cerr;
extern MasterStream clog;

}
#endif
