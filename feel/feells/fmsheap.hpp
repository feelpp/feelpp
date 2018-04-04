/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#ifndef _FMSHEAP_HPP
#define _FMSHEAP_HPP 1

#include <functional>

namespace Feel
{
namespace details
{

struct none_type
{};

/* 
 *  FmsHeap is a heap structure, used to retrieve the min of abs(phi) when performing fast-marching.
 */
template< typename T, typename DataType = none_type >
class FmsHeap;

// Generic heap (with data)
template< typename T, typename DataType >
class FmsHeap
{
public:
    typedef T value_type;
    typedef DataType data_type;

    // typically : pair < phi_value,  index >
    typedef std::pair<value_type, size_type> heap_entry_type;

    struct heap_entry_data_type
        : public std::pair<heap_entry_type, data_type>
    {
        typedef std::pair<heap_entry_type, data_type> super_type;
        // Inherit constructor
        using typename super_type::pair;
        heap_entry_data_type(heap_entry_type const& entry, data_type const& v)
            : super_type(entry, v)
        {}
        heap_entry_data_type() : super_type() {}

        heap_entry_type & heap_entry() { return this->first; }
        heap_entry_type const& heap_entry() const { return this->first; }

        value_type & value() { return this->first.first; }
        value_type const& value() const { return this->first.first; }

        size_type & index() { return this->first.second; }
        size_type const& index() const { return this->first.second; }

        data_type & data() { return this->second; }
        data_type const& data() const { return this->second; }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->first;
            ar & this->second;
        }
    };

    //--------------------------------------------------------------------//
    // Constructor
    FmsHeap() : maxValue( std::numeric_limits<value_type>::max() ) {}

    void push( heap_entry_type in, data_type const& data = data_type() )
    {
        M_heap.push_back( heap_entry_data_type(in, data) );
        push_heap( M_heap.begin(), M_heap.end(), cmp );
    }

    void change( heap_entry_type in, data_type const& data )
    {
        for ( auto it = M_heap.begin(); it != M_heap.end(); ++it )
        {
            // if the entry at index exists,
            // this phi_entry take the min value between |phi_entry| and |phi_min|
            if ( it->index() == in.second )
            {
                if ( cmp_entry( it->heap_entry(), in ) )
                {
                    *it = heap_entry_data_type(in, data);
                    make_heap( M_heap.begin(), M_heap.end(), cmp );
                }
                return;
            }
        }
        // if not found, push
        this->push( in, data );
    }
    void change( heap_entry_type in )
    {
        for ( auto it = M_heap.begin(); it != M_heap.end(); ++it )
        {
            // if the entry at index exists,
            // this phi_entry take the min value between |phi_entry| and |phi_min|
            if ( it->first.second == in.second )
            {
                if ( cmp_entry( it->first, in ) )
                {
                    *it = heap_entry_data_type(in, it->data());
                    make_heap( M_heap.begin(), M_heap.end(), cmp );
                }
                return;
            }
        }
        // if not found, push
        this->push( in, data_type() );
    }

    heap_entry_data_type pop()
    {
        // assert M_heap.size() > 0
        heap_entry_data_type out = *(M_heap.begin());
        pop_heap( M_heap.begin(), M_heap.end(), cmp );
        M_heap.pop_back();
        return out;
    }


    /* returns the front element, which is the one having the smallest phi
       if the heap is empty, return an element having a huge value of phi, so that, compared to other heap entries it will never be accepted */
    heap_entry_data_type front()
    {
        if ( M_heap.empty() )
            return heap_entry_data_type(std::make_pair(maxValue, 0 ), data_type());
        return M_heap.front();
    }
    heap_entry_type frontEntry()
    {
        return this->front().heap_entry();
    }
    data_type frontData()
    {
        return this->front().data();
    }

    bool checkExistingEntry(size_type index)
    {
        for (auto it = M_heap.begin(); it != M_heap.end(); ++it )
                if (it->index() == index)
                    return true;
        return false;
    }

    /* remove from the heap the entry having the specific value index*/
    bool removeFromHeap(size_type index)
    {
        bool removed = false;
        for (auto it = M_heap.begin(); it != M_heap.end(); ++it )
            if (it->first.second == index)
            {
                M_heap.erase(it);
                removed = true;
                break;
            }
        return removed;
    }

    static heap_entry_type min(heap_entry_type a, heap_entry_type b)
    {
        return std::abs(a.first < b.first) ? a : b;
    }

    typename std::vector<heap_entry_type>::iterator begin() { return M_heap.begin(); }
    typename std::vector<heap_entry_type>::iterator end() { return M_heap.end(); }
    size_type size() const { return M_heap.size(); }

    void clear() { M_heap.clear(); }

    value_type valueAtIndex( size_type index )
    {
        for (auto const& heapEntry : M_heap)
            if (heapEntry.index() == index )
                return heapEntry.value();

        CHECK( false ) << "index: "<<index<<" does not exists in the heap\n";
        return 0;
    }

    data_type const& dataAtIndex( size_type index ) const
    {
        for (auto const& heapEntry : M_heap)
            if (heapEntry.index() == index )
                return heapEntry.data();

        CHECK( false ) << "index: "<<index<<" does not exists in the heap\n";
        return data_type();
    }

    data_type & dataAtIndex( size_type index )
    {
        for (auto & heapEntry : M_heap)
            if (heapEntry.index() == index )
                return heapEntry.data();

        CHECK( false ) << "index: "<<index<<" does not exists in the heap\n";
        //return data_type();
    }

private:
    const value_type maxValue;

    static bool cmp_entry( heap_entry_type a, heap_entry_type b )
    {
        // returns true if |phi_a| > |phi_b|
        value_type aa = a.first;
        value_type bb = b.first;
        return std::abs(aa) > std::abs(bb);
    }
    static bool cmp( heap_entry_data_type a, heap_entry_data_type b )
    {
        return cmp_entry(a.first, b.first);
    }

    std::vector<heap_entry_data_type> M_heap;
};

// Specilization for no-data case
template<typename T>
class FmsHeap<T, none_type>
{
public:
    typedef T value_type;
    typedef none_type data_type;

    //         pair < phi_value,  index >
    typedef std::pair<value_type, size_type> heap_entry_type;
    struct heap_entry_data_type
        : public heap_entry_type
    {
        typedef heap_entry_type super_type;
        // Inherit constructor
        using typename super_type::pair;
        heap_entry_data_type(heap_entry_type const& entry, data_type const& = data_type() )
            : super_type(entry)
        {}
        heap_entry_data_type() : super_type() {}

        heap_entry_type & heap_entry() { return *this; }
        heap_entry_type const& heap_entry() const { return *this; }

        value_type & value() { return this->first; }
        value_type const& value() const { return this->first; }

        size_type & index() { return this->second; }
        size_type const& index() const { return this->second; }

        //data_type & data() { return this->second; }
        data_type data() const { return data_type(); }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<super_type>(*this);
        }
    };

    //--------------------------------------------------------------------//
    // Constructor
    FmsHeap() : maxValue( std::numeric_limits<value_type>::max() ) {}

    void push( heap_entry_type in, data_type const& data = data_type() )
    {
        M_heap.push_back( in );
        push_heap( M_heap.begin(), M_heap.end(), cmp );
    }

    void change( heap_entry_type in, data_type const& data )
    {
        this->change(in);
    }
    void change( heap_entry_type in )
    {
        for ( auto it = M_heap.begin(); it != M_heap.end(); ++it )
        {
            // if the entry at index exists,
            // this phi_entry take the min value between |phi_entry| and |phi_min|
            if ( it->index() == in.second )
            {
                if ( cmp_entry( it->heap_entry(), in ) )
                {
                    *it = in;
                    make_heap( M_heap.begin(), M_heap.end(), cmp );
                }
                return;
            }
        }
        // if not found, push
        this->push( in );
    }

    heap_entry_data_type pop()
    {
        // assert M_heap.size() > 0
        heap_entry_data_type out = *(M_heap.begin());
        pop_heap( M_heap.begin(), M_heap.end(), cmp );
        M_heap.pop_back();
        return out;
    }


    /* returns the front element, which is the one having the smallest phi
       if the heap is empty, return an element having a huge value of phi, so that, compared to other heap entries it will never be accepted */
    heap_entry_data_type front()
    {
        if ( M_heap.empty() )
            return heap_entry_data_type(std::make_pair(maxValue, 0 ));
        return M_heap.front();
    }
    heap_entry_type frontEntry()
    {
        return this->front().heap_entry();
    }
    data_type frontData()
    {
        return this->front().data();
    }

    bool checkExistingEntry(size_type index)
    {
        for (auto it = M_heap.begin(); it != M_heap.end(); ++it )
            if (it->index() == index)
                return true;
        return false;
    }

    /* remove from the heap the entry having the specific value index*/
    bool removeFromHeap(size_type index)
    {
        bool removed = false;
        for (auto it = M_heap.begin(); it != M_heap.end(); ++it )
            if (it->index() == index)
            {
                M_heap.erase(it);
                removed = true;
                break;
            }
        return removed;
    }

    static heap_entry_type min(heap_entry_type a, heap_entry_type b)
    {
        return std::abs(a.first < b.first) ? a : b;
    }

    typename std::vector<heap_entry_type>::iterator begin() { return M_heap.begin(); }
    typename std::vector<heap_entry_type>::iterator end() { return M_heap.end(); }
    size_type size() const { return M_heap.size(); }

    void clear() { M_heap.clear(); }

    value_type valueAtIndex( size_type index )
    {
        for (auto const& heapEntry : M_heap)
            if (heapEntry.index() == index )
                return heapEntry.value();

        CHECK( false ) << "index: "<<index<<" does not exists in the heap\n";
        return 0;
    }

    data_type dataAtIndex( size_type index ) const
    {
        return data_type();
    }

    //data_type & dataAtIndex( size_type index )
    //{
        //for (auto & heapEntry : M_heap)
            //if (heapEntry.index() == index )
                //return heapEntry.data();

        //CHECK( false ) << "index: "<<index<<" does not exists in the heap\n";
        ////return data_type();
    //}

private:
    const value_type maxValue;

    static bool cmp_entry( heap_entry_type a, heap_entry_type b )
    {
        // returns true if |phi_a| > |phi_b|
        value_type aa = a.first;
        value_type bb = b.first;
        return std::abs(aa) > std::abs(bb);
    }
    static bool cmp( heap_entry_data_type a, heap_entry_data_type b )
    {
        return cmp_entry(a, b);
    }

    std::vector<heap_entry_data_type> M_heap;
};

} // namespace details

} // namespace Feel

#endif /* _FMSHEAP_HPP */
