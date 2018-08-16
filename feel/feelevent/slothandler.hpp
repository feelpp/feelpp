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

#ifndef FEELPP_SLOTHANDLER_HPP
#define FEELPP_SLOTHANDLER_HPP 1

#include <iostream>
#include <map>
#include <functional>
#include <boost/any.hpp>

namespace Feel
{
namespace Event
{

class SlotHandler
{
public:

    //! Feel++ slot type.
    template< typename T >
        using func_type = std::function<T>;
    using slot_map_type = std::map< std::string, boost::any >;
    template< typename SlotType >
        using slot_map_entry_type = std::pair< std::string,
                                               func_type<SlotType> >;

    //| Allocators
    //! @{

    //! Create a new slot.
    //! \tparam Slotype Slot functor type to apply to a given signal.
    //! \param name Name of the slot
    //! \param slot A functor called by a signal. (default empty)
    //! By default, the slot is empty and has to be set via slotSet in that case.
    //! \return A pointer to the stored slot (Return a nullptr if no functor is
    //!  passed as argument).
    //! \see slotSet
//    template< typename SlotType >
//        const auto& slotNew( const std::string& name )
//    {
//        return slotStore<boost::any>( name, slot_any_type() );
//    };

    //! Remove a slot.
    void slotDelete( const std::string& name )
    {
        M_slots.erase( name );
    };

    //! Create a new slot.
    //! \tparam Slotype Slot functor type to apply to a given signal.
    //! \param name Name of the slot
    //! \param slot A functor called by a signal. (default empty)
    //! By default, the slot is empty and has to be set via slotSet in that case.
    //! \return A pointer to the stored slot (Return a nullptr if no functor is
    //!  passed as argument).
    //! \see slotSet
    template< typename SlotType >
        const auto& slotNew( const std::string& name,
                             const func_type<SlotType>& slot )
    {
        return slotStore( name, slot );
    };

    template< typename SlotType >
        static const auto& slotStaticNew( const std::string& name,
                             const func_type<SlotType>& slot )
    {
        return slotStaticStore( name, slot );
    };

    //! @}

    //! Setters
    //! @{
    //! @}

    //! Getters
    //! @{

    //! Get a slot by name (any).
    //! Note: This function return a type "any". The user has to perform a
    //! dynamic cast afterwards to the real slot type.
    const boost::any& slot( const std::string& name )
    {
        return M_slots[name];
    };

    //! Get a static slot by name (any).
    //! Note: This function return a type "any". The user has to perform a
    //! dynamic cast afterwards to the real slot type.
    static const boost::any& slotStatic( const std::string& name )
    {
        return S_slots[name];
    };

    //! Get a slot by name (slot type).
    //! Note: This function return the real slot type passed by template
    //! argument.
    template< typename SlotType >
    decltype(auto)
    slot( const std::string& name )
    {
        using RetType = const func_type<SlotType>&;
        return boost::any_cast< RetType >( M_slots[name] );
    };

    //! Get a slot by name (slot type).
    //! Note: This function return the real slot type passed by template
    //! argument.
    template< typename SlotType >
    static decltype(auto)
    slotStatic( const std::string& name )
    {
        using RetType = const func_type<SlotType>&;
        return boost::any_cast< RetType >( S_slots[name] );
    };

    // Get a map of all signals.
    const slot_map_type& slots() const
    {
        return M_slots;
    };

    // Get a static map of all signals.
    static const slot_map_type& slotsStatic()
    {
        return S_slots;
    };

    // @}

    void slotShow() const
    {
        std::cout << "SLOTS:" << std::endl;
        std::cout << std::string(40,'-') << std::endl;
        for( const auto& l : M_slots )
            std::cout << "* " << l.first << std::endl;
        std::cout << std::string(40,'-') << std::endl;
    }

    void slotStaticShow() const
    {
        std::cout << "STATIC SLOTS:" << std::endl;
        std::cout << std::string(40,'-') << std::endl;
        for( const auto& l : S_slots )
            std::cout << "* " << l.first << std::endl;
        std::cout << std::string(40,'-') << std::endl;
    }

    //! Operators
    //! @{

    // Compare operators.
    friend bool operator< (const SlotHandler& lhs, const SlotHandler& rhs )
    {
        return &lhs < &rhs;
    }

    //! @}

private:

    //! Store a slot as a pointer in the map.
    template< typename SlotType >
    decltype(auto)
    slotStore( const std::string& name,
               const func_type<SlotType>& s )
    {
        using RetType = const func_type<SlotType> &;
        using MapEntry = slot_map_entry_type< SlotType >;
        M_slots.insert( MapEntry( name, s ) );
        return boost::any_cast< RetType >( M_slots[name] );
    }

    //! Store a slot as a pointer in the map.
    template< typename SlotType >
    static decltype(auto)
    slotStaticStore( const std::string& name,
               const func_type<SlotType>& s )
    {
        using RetType = const func_type<SlotType> &;
        using MapEntry = slot_map_entry_type< SlotType >;
        S_slots.insert( MapEntry( name, s ) );
        return boost::any_cast< RetType >( S_slots[name] );
    }

private:
    //! Map containing signal name and object.
    slot_map_type M_slots;
    static slot_map_type S_slots;
};

} // Event namespace
} // Feel namespace.

#endif // FEELPP_SLOTHANDLER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
