/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-23

  Copyright (C) 2007,2009 Université de Grenoble 1

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file timermap.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-23
 */
#ifndef __TimerMap_H
#define __TimerMap_H 1

#include <boost/timer.hpp>

namespace Feel
{
struct TimerPair : public std::pair<boost::timer, double>
{
public:

    double elapsed()
    {
        second = first.elapsed();
        return second;
    }
    void accumulate()
    {
        second += first.elapsed();
        first.restart();
    }
    void restart()
    {
        first.restart();
    }
    void reset()
    {
        first.restart();
        second = 0;
    }
};
/**
 * \class TimerMap
 * \brief timers map
 *
 * @author Christophe Prud'homme
 * @see
 */
class TimerMap : public std::map<std::string, TimerPair>
{
    typedef std::map<std::string, TimerPair> super;
public:



    /** @name Typedefs
     */
    //@{

    typedef super::iterator iterator;
    typedef super::const_iterator const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    TimerMap() {}
    TimerMap( TimerMap const & tm ) : super( tm ) {}
    ~TimerMap() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    template<typename PrefixType>
    void report( PrefixType const& prefix )
    {
        const_iterator it;

        for ( it=this->begin(); it!=this->end(); ++it )
        {
            VLOG(1) << prefix << " "  << it->first << ": " << it->second.second << "\n";
        }
    }

    //@}



protected:

private:

};
} // Feel
#endif /* __TimerMap_H */
