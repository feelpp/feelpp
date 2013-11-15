/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-10-23

  Copyright (C) 2005,2006 EPFL

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
   \file fmsheap.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-10-23
 */
#ifndef __FMS_Heap_H
#define __FMS_Heap_H 1

namespace Feel
{
namespace details
{

template<typename T>
class FmsHeap
{

public:

    typedef T value_type;
    typedef std::pair<value_type, uint16_type> heap_entry_type;

    void push( heap_entry_type in )
    {
        M_heap.push_back( in );
        push_heap( M_heap.begin(), M_heap.end(), farther );
    }

    void change( heap_entry_type in )
    {
        for ( heapvect_iter it=M_heap.begin(); it != M_heap.end(); ++it )
        {
            if ( it->second == in.second )
            {
                if ( farther( *it, in ) )
                {
                    *it = in;
                    make_heap( M_heap.begin(),
                               M_heap.end(),
                               farther );
                }

                return;
            }
        }

        // if not found, push
        push( in );
    }

    heap_entry_type pop()
    {
        // assert M_heap.size() > 0
        heap_entry_type out = *( M_heap.begin() );
        pop_heap( M_heap.begin(), M_heap.end(), farther );
        M_heap.pop_back();
        return out;
    }

    size_type size() const
    {
        return M_heap.size();
    }

private:

    typedef std::vector<heap_entry_type> heapvect_type;
    typedef typename heapvect_type::iterator heapvect_iter;
    typedef typename heapvect_type::const_iterator heapvect_constiter;

    static bool farther( heap_entry_type a, heap_entry_type b )
    {
        value_type aa = a.first;
        aa = aa < 0.0 ? -aa : aa;
        value_type bb = b.first;
        bb = bb < 0.0 ? -bb : bb;
        return aa > bb;
    }

    std::vector<heap_entry_type> M_heap;
};

} // namespace details

} // namespace Feel

#endif /* __FMS_Heap_H */
