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

#ifndef FEELPP_HWSYSBASE_HPP
#define FEELPP_HWSYSBASE_HPP 1

#include <feel/feelobserver/journalwatcher.hpp>

namespace Feel
{

static uint16_type HWSYS_INSTANCE_NUMBER = 0;

namespace Sys
{

namespace pt =  boost::property_tree;

//! Hardware Sytem base class.
//! This class gather information from the operating system such
//! as memory, number of processors, distribution ...
//! External libraries should have a specific header inheriting from
//! this base class.
class HwSysBase
    : public Observer::JournalWatcher
{
public:

    //! Default constructor.
    HwSysBase()
        : M_backend( "hwsys" )
    {
        instanceName( "hwsys-" + std::to_string(HWSYS_INSTANCE_NUMBER) );
        HWSYS_INSTANCE_NUMBER++;
        // TODO fix/remove this line. This should not be necessary!
        // There's must be a problem with the inheritance here.
        this->journalDisconnect();
        VLOG(2) << "[HwSys] constructor instance number " << M_instance_number << "\n" ;
    }

    //! Copy constructor.
    HwSysBase( HwSysBase& ) = default;

    //! Move constructor.
    HwSysBase( HwSysBase&& ) = default;

    //! Assignement constructor.
    HwSysBase& operator=( HwSysBase& ) = default;
    
    //! Assignement-move constructor.
    HwSysBase& operator=( HwSysBase&& ) = default;

    //! Destructor
    virtual ~HwSysBase() override = default;

    //! Accessors
    //! @{

    //! Get the name of the current mesh instance.
    //! @return the name
    std::string instanceName() const { return M_instance_name; }

    const std::string& instanceName( std::string s ) 
    {
        M_instance_name=s;
        return M_instance_name;
    }

#if 0
    const std::string& OSName( std::string const& s ) { return M_os_name; }
    const std::string& OSRelease( std::string const& s ) { return M_os_release; }
    const std::string& OSVersion( std::string const& s ) { return M_os_version; }
    const std::string& OSPlatform( std::string const& s ) { return M_os_platform; }
    const std::string& HostName( std::string const& s ) { return M_host_name; }
    const std::string& DomainName( std::string const& s ) { return M_domain_name; }

    const std::string& HostTotalMemory( std::string const& s ) { return M_host_total_mem; }
    const std::string& HostTotalMemoryAvailable( std::string const& s ) { return M_host_total_mem_avail; }
    const std::string& ProcCacheSize( std::string const& s ) { return M_proc_cache_size; }
#endif

    // @}
    
private: 

    //! Private Methods
    //! @{
    
    //! Journal watcher notification.
    const pt::ptree journalNotify() const override
    {
        pt::ptree p;
        const auto ccp = std::to_string( Observer::JournalManager::journalCurrentCheckpoint() );
		if( not M_host_name.empty() )
		{
            std::string prefix = "hardware." + M_host_name;

		    p.put( prefix + ".backend", M_backend );
			p.put( prefix + ".os.name", M_os_name );
			p.put( prefix + ".is_linux", M_os_is_linux );
			p.put( prefix + ".is_apple", M_os_is_apple );
			p.put( prefix + ".is_windows", M_os_is_windows );
			p.put( prefix + ".os.release", M_os_release );
			p.put( prefix + ".os.version", M_os_version );
			p.put( prefix + ".os.platform", M_os_platform );
			p.put( prefix + ".hostname", M_host_name );
			p.put( prefix + ".domain_fqname", M_domain_name );
			p.put( prefix + ".is_64bits", M_proc_is64bits );
			p.put( prefix + ".proc.vendor_name", M_proc_vendor_name );
			p.put( prefix + ".proc.vendor_id", M_proc_vendor_id );
			p.put( prefix + ".proc.type_id", M_proc_type_id );
			p.put( prefix + ".proc.family_id", M_proc_family_id );
			p.put( prefix + ".proc.model_id", M_proc_model_id );
			p.put( prefix + ".proc.extended_name", M_proc_extended_name );
			p.put( prefix + ".proc.stepping_code", M_proc_stepping_code );
			p.put( prefix + ".proc.serial_number", M_proc_serial_number );
			p.put( prefix + ".proc.cache_size", M_proc_cache_size );
			p.put( prefix + ".proc.logical_per_physical", M_proc_logical_per_physical );
			p.put( prefix + ".proc.clock_frequency", M_proc_clock_frequency );
			p.put( prefix + ".proc.logical_cpu_number", M_proc_logical_cpu_number );
			p.put( prefix + ".proc.physical_cpu_number", M_proc_physical_cpu_number );
			p.put( prefix + ".proc.cpu_id_support", M_proc_cpu_id_support );
			p.put( prefix + ".proc.apic_id", M_proc_apic_id );
			p.put( prefix + ".mem.total.virtual", M_mem_virtual_total );
			p.put( prefix + ".mem.total.physical", M_mem_physical_total );
			p.put( prefix + ".mem.total.host", M_mem_host_total );
			p.put( prefix + ".mem.available.checkpoint-" + ccp + ".virtual", M_mem_virtual_avail );
			p.put( prefix + ".mem.available.checkpoint-" + ccp + ".physical", M_mem_physical_avail );
			p.put( prefix + ".mem.available.checkpoint-" + ccp + ".host", M_mem_host_avail );
			p.put( prefix + ".mem.available.checkpoint-" + ccp + ".proc", M_mem_proc_avail );
			p.put( prefix + ".mem.used.checkpoint-" + ccp + ".host", M_mem_host_used );
			p.put( prefix + ".mem.used.checkpoint-" + ccp + ".proc", M_mem_proc_used );
			p.put( prefix + ".mem.load_average", M_load_avg );
		}

        return p;
    }
    
    //! @}

protected:

    std::string M_instance_name;
    std::string M_instance_number;
    std::string M_backend;

    //! Operating System
    //! @{
    std::string M_os_name;
    std::string M_os_is_linux;
    std::string M_os_is_apple;
    std::string M_os_is_windows;
    std::string M_os_release;
    std::string M_os_version;
    std::string M_os_platform;

    std::string M_host_name;
    std::string M_domain_name;
    //! @}

    // Hardware
    // @{
    std::string M_proc_is64bits;
    std::string M_proc_vendor_name;
    std::string M_proc_vendor_id;
    std::string M_proc_type_id;
    std::string M_proc_family_id;
    std::string M_proc_model_id;
    std::string M_proc_extended_name;
    std::string M_proc_stepping_code;
    std::string M_proc_serial_number;
    std::string M_proc_cache_size;
    std::string M_proc_logical_per_physical;
    std::string M_proc_clock_frequency;
    std::string M_proc_logical_cpu_number;
    std::string M_proc_physical_cpu_number;
    std::string M_proc_cpu_id_support;
    std::string M_proc_apic_id;

	std::string M_mem_virtual_total;
	std::string M_mem_virtual_avail;
	std::string M_mem_physical_total;
	std::string M_mem_physical_avail;
	std::string M_mem_host_total;
	std::string M_mem_host_avail;
	std::string M_mem_proc_avail;
	std::string M_mem_host_used;
	std::string M_mem_proc_used;
	std::string M_load_avg;

    //! @}
};


} // Sys
} // Feel


#endif // FEELPP_HWSYSBASE_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
