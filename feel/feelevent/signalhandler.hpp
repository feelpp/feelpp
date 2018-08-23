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

#ifndef FEELPP_SIGNALHANDLER_HPP
#define FEELPP_SIGNALHANDLER_HPP 1

#include <boost/signals2.hpp>
#include <boost/any.hpp>
#include <feel/feelcore/feelassert.hpp>
#include <feel/feelevent/slothandler.hpp>

namespace Feel
{
namespace Event
{

//
enum SigConnectFlags{
    DEFAULT,
    SIG_STATIC,
    SLOT_STATIC
};

inline SigConnectFlags operator| (SigConnectFlags a, SigConnectFlags b)
{
    return static_cast<SigConnectFlags>( static_cast<int>(a) | static_cast<int>(b) );
}

//! Enable signals & slots to inherited objects.
//! Add new functionnalities to connect objects between each others.
//! Class that inherite from SignalHandler will have to possibility to create
//! new signal.
//!
//! :Note: Function template is a design choice. We want to simplify event
//! creation in userland and also minimize polymorphism complexity.
//!
//! :Info: This class is used to generate benchmark database.
class SignalHandler
{
public:

    //! Typedefs
    //! @{

    //! Feel++ signal type.
    //! \tparam  Args Signal arguments
    //! The first argument is the Slot type, second is the event to apply to all
    //! notification.
    template< typename... Args >
        using sig_type = boost::signals2::signal< Args... >;
    template< typename... Args >
        using sig_shared_ptr_type = std::shared_ptr< sig_type< Args... > >;
    template< typename... Args >
        using sig_weak_ptr_type = std::shared_ptr< sig_type< Args... > >;
    //| Map types used to store the signal.
    using sig_map_type = std::map< std::string, boost::any >;
    template< typename... Args >
        using sig_map_entry_type = std::pair< std::string,
                                              sig_shared_ptr_type<Args...> >;
    using slothdlr_base_type = SlotHandler;
    using connection_type = boost::signals2::connection;
    using link_type = std::tuple<std::string, std::reference_wrapper<slothdlr_base_type const>, std::string>;
    using link_map_type = std::map< link_type, std::vector<connection_type> >;

    //! @}

    // Allocators
    // @{

    //! Create a new signal.
    //! \tparam  Args optional signal arguments (see sig_type).
    template< typename... Args >
    const sig_shared_ptr_type< Args... >&
    signalNew( const std::string& name )
    {
        using SigType = sig_type<Args...>;
        const auto& sig = std::make_shared< SigType >( SigType() );
        return signalStore<Args...>( name, sig );
    }

    //! Create a new static signal.
    //! \tparam  Args optional signal arguments (see sig_type).
    template< typename... Args >
    static const sig_shared_ptr_type< Args... >&
    signalStaticNew( const std::string& name )
    {
        using SigType = sig_type<Args...>;
        const auto& sig = std::make_shared< SigType >( SigType() );
        return signalStaticStore<Args...>( name, sig );
    }

    //! Delete a signal.
    void signalDelete( const std::string& name )
    {
        M_sigs.erase( name );
    }

    //! Delete a signal.
    static void signalStaticDelete( const std::string& name )
    {
        S_sigs.erase( name );
    }

    // @}

    //! Getters
    //! @{

    //! Get a signal by name.
    //! Use dynamic cast to get real signal type.
    //! Otherwise, prefer using createSignal return type directly.
    //! \return a signal as any type. A dynamic_cast has to be performed
    //! afterward to retrieve real signal type!
    const boost::any& signal( const std::string& name )
    {
        return M_sigs[name];
    }

    //! Get a signal by name.
    //! Otherwise, prefer using createSignal return type directly.
    //! \return the real signal type. The slot type has
    //! to be passed as a template argument.
    template< typename SlotType, typename... Args >
    decltype(auto)
    signal( const std::string& name )
    {
        using RetType = const SignalHandler::sig_shared_ptr_type< SlotType, Args... >&;
        const boost::any& asig = M_sigs.at(name);
        RetType sig = boost::any_cast< RetType >( asig );
        CHECK( sig != nullptr ) << "Bad signal allocation 'signalHandler'";
        return sig;
    }

    //! Get a static signal by name.
    //! Otherwise, prefer using createSignal return type directly.
    //! \return the real signal type. The slot type has
    //! to be passed as a template argument.
    template< typename SlotType, typename... Args >
    static decltype(auto)
    signalStatic( const std::string& name )
    {
        using RetType = const SignalHandler::sig_shared_ptr_type< SlotType, Args... >&;
        const boost::any& asig = S_sigs.at(name);
        RetType sig = boost::any_cast< RetType >( asig );
        CHECK( sig != nullptr ) << "Bad signal allocation 'SignalHandler'";
        return  sig;
    }

    //! Get a signal by name.
    //! Use dynamic cast to get real signal type.
    //! Otherwise, prefer using createSignal return type directly.
    const sig_map_type& signals() const
    {
        return M_sigs;
    }

    //! Get a static signal by name.
    //! Use dynamic cast to get real signal type.
    //! Otherwise, prefer using createSignal return type directly.
    static const sig_map_type& signalsStatic()
    {
        return S_sigs;
    }

    // @}

    //! Connect a slot to a signal.
    //! \tparam SlotType Slot type.
    //! \tparam SignalArgs Signal optional argument types.
    //! \tparam SlotHdlr Slot handler type.
    //! \param signame Signal name.
    //! \param slothdlr Slot handler.
    //! \param slotname Slot name.
    //! The first argument is the slot type. The signal can only connect
    //! to a slot that have the same type!
    template< typename SlotType, typename... SignalArgs, typename SlotHdlr >
    connection_type signalConnect( const std::string& signame,
                                   SlotHdlr& slothdlr,
                                   const std::string& slotname,
                                   SigConnectFlags flags = DEFAULT )
    {
        const auto& sig = signal<SlotType, SignalArgs...>( signame );
        const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
        if( flags == SLOT_STATIC )
        {
            const auto& slo = slothdlr.template slotStatic<SlotType>( slotname );
            const auto& c = sig->connect( slo );
            S_links[link].push_back( c );
        }
        else
        {
            const auto& slo = slothdlr.template slot<SlotType>( slotname );
            const auto& c = sig->connect( slo );
            S_links[link].push_back( c );
        }
    }

    //! Connect a slot to a static signal.
    //! \tparam SlotType Slot type.
    //! \tparam SignalArgs Signal optional argument types.
    //! \tparam SlotHdlr Slot handler type.
    //! \param signame Signal name.
    //! \param slothdlr Slot handler.
    //! \param slotname Slot name.
    //! The first argument is the slot type. The signal can only connect
    //! to a slot that have the same type!
    template< typename SlotType, typename... SignalArgs, typename SlotHdlr >
    static const connection_type
    signalStaticConnect( const std::string& signame,
                         SlotHdlr& slothdlr,
                         const std::string& slotname,
                         SigConnectFlags flags = DEFAULT )
    {
       const auto& sig = signalStatic<SlotType, SignalArgs...>( signame );
       const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
       if( flags == SLOT_STATIC )
       {
           const auto& slo = slothdlr.template slotStatic<SlotType>( slotname );
           const auto& c = sig->connect( slo );
           S_links[link].push_back( c );
       }
       else
       {
           const auto& slo = slothdlr.template slot<SlotType>( slotname );
           const auto& c = sig->connect( slo );
           S_links[link].push_back( c );
       }
       return S_links[link].back();
    }

    //! Disconnect a slot from the signal.
    //! \tparam SlotType Slot type.
    //! \tparam SignalArgs Signal optional argument types.
    //! \tparam SlotHdlr Slot handler type.
    //! \param signame Signal name.
    //! \param slothdlr Slot handler.
    //! \param slotname Slot name.
    //! The first argument is the slot type. The signal can only connect
    //! to a slot that have the same type!
    template< typename SlotType, typename... SignalArgs, typename SlotHdlr >
    void signalDisconnect( const std::string& signame,
                        SlotHdlr& slothdlr,
                        const std::string& slotname,
                        SigConnectFlags flags = DEFAULT )
    {
        const auto& sig = signal<SlotType, SignalArgs...>( signame );
        if( flags == SLOT_STATIC )
        {
            const auto& slo = slothdlr.template slot<SlotType>( slotname );
            const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
            const auto& v = S_links[link];
            for( const auto& c : v )
                c.disconnect();
        }
        else
        {
            const auto& slo = slothdlr.template slot<SlotType>( slotname );
            const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
            const auto& v = S_links[link];
            for( const auto& c : v )
                c.disconnect();
        }
    }

    //! Disconnect a slot from a static signal.
    //! \tparam SlotType Slot type.
    //! \tparam SignalArgs Signal optional argument types.
    //! \tparam SlotHdlr Slot handler type.
    //! \param signame Signal name.
    //! \param slothdlr Slot handler.
    //! \param slotname Slot name.
    //! The first argument is the slot type. The signal can only connect
    //! to a slot that have the same type!
    //! Note: All connections are closed!
    template< typename SlotType, typename... SignalArgs, typename SlotHdlr >
    static void
    signalStaticDisconnect( const std::string& signame,
                         SlotHdlr& slothdlr,
                         const std::string& slotname,
                         SigConnectFlags flags = DEFAULT )
    {
       const auto& sig = signalStatic<SlotType, SignalArgs...>( signame );
       if( flags == SLOT_STATIC )
       {
           const auto& slo = slothdlr.template slotStatic<SlotType>( slotname );
           const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
           const auto& v = S_links[link];
           for( const auto& c : v )
               c.disconnect();
       }
       else
       {
           const auto& slo = slothdlr.template slot<SlotType>( slotname );
           const link_type& link = std::make_tuple( signame, std::cref(slothdlr), slotname );
           const auto& v = S_links[link];
           for( const auto& c : v )
               c.disconnect();
       }
    }

    //! Print a list of signal.
    //! Display stored signal from the map.
    void signalShow() const
    {
        std::cout << "SIGNALS:" << std::endl;
        std::cout << std::string(40,'-') << std::endl;
        for( const auto& l : M_sigs )
            std::cout << "* " << l.first << std::endl;
        std::cout << std::string(40,'-') << std::endl;
    }

    //! Print a list of static signal.
    //! Display stored signal from the static map.
    void signalStaticShow() const
    {
        std::cout << "STATIC SIGNALS:" << std::endl;
        std::cout << std::string(40,'-') << std::endl;
        for( const auto& l : S_sigs )
            std::cout << "* " << l.first << std::endl;
        std::cout << std::string(40,'-') << std::endl;
    }

    //! Fetch notifications coming from slots connected to a given signal.
    template< typename... SignalArgs >
    decltype(auto)
    signalPull( std::string name )
    {
       const auto& sig = signal<SignalArgs...>( name );
       return (*sig)();
    }

    //! Fetch notifications coming from slots connected to a given signal.
    template< typename... SignalArgs >
    static decltype(auto)
    signalStaticPull( std::string name )
    {
       const auto& sig = signalStatic<SignalArgs...>( name );
       return (*sig)();
    }

private:

    //! Store a signal into the map and return pointer
    //! to this stored signal.
    template< typename... Args >
    decltype(auto)
    signalStore( std::string name,
                 const sig_shared_ptr_type<Args...>& sig )
    {
        // Check that the sig does not exist in the static map.
        using RetType = const SignalHandler::sig_shared_ptr_type< Args... >&;
        using MapEntry = SignalHandler::sig_map_entry_type<Args...>;
        M_sigs.insert( MapEntry( name, sig ) );
        return boost::any_cast< RetType >( M_sigs[name] );
    }

    //! Store a signal into the static map and return a pointer
    //! to this stored signal.
    template< typename... Args >
    static decltype(auto)
    signalStaticStore( std::string name,
                       const sig_shared_ptr_type<Args...>& sig )
    {
        // Check that the sig does not exist in the non static map.
        using RetType = const SignalHandler::sig_shared_ptr_type< Args... >&;
        using MapEntry = SignalHandler::sig_map_entry_type<Args...>;
        S_sigs.insert( MapEntry( name, sig ) );
        return boost::any_cast< RetType >( S_sigs[name] );
    }

private:

    //! Map containing signal..
    sig_map_type M_sigs;
    //! Static Map containing signal
    static sig_map_type S_sigs;
    //! Static map gathering connections per link.
    static link_map_type S_links;
};

} // Event namespace
} // Feel namespace.

#endif // FEELPP_SIGNALHANDLER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
