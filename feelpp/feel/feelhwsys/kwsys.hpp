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

#ifndef FEELPP_KWSYS_HPP
#define FEELPP_KWSYS_HPP 1

#include <feel/feelhwsys/hwsysbase.hpp>

// contrib headers.
#include <feelpp_kwsys/SystemInformation.hxx>

namespace Feel
{
namespace Sys
{

//! Wrapper class for Kitware system library (kwsys).
class KWSys
    : public HwSysBase
{
public:
    KWSys()
    {
        using namespace feelpp_kwsys;
        SystemInformation info;
        
        // Use KWSys.
        M_backend = "kwsys";

        info.RunCPUCheck();
        info.RunOSCheck();
        info.RunMemoryCheck();

        M_os_name = info.GetOSName();
        M_os_is_linux = std::to_string( info.GetOSIsLinux() );
        M_os_is_apple = std::to_string( info.GetOSIsApple() );
        M_os_is_windows = std::to_string( info.GetOSIsWindows() );
        M_os_release = info.GetOSRelease();
        M_os_version = info.GetOSVersion();
        M_os_platform = info.GetOSPlatform();
        M_host_name = info.GetHostname();
        M_domain_name = info.GetFullyQualifiedDomainName();

        M_proc_is64bits = std::to_string( info.Is64Bits() );
        M_proc_vendor_name = info.GetVendorString();
        M_proc_vendor_id = info.GetVendorID();
        M_proc_type_id = info.GetTypeID();
        M_proc_family_id = info.GetFamilyID();
        M_proc_model_id = info.GetModelID();
        M_proc_extended_name = info.GetExtendedProcessorName();
        M_proc_stepping_code = info.GetSteppingCode();
        M_proc_serial_number = info.GetProcessorSerialNumber();
        
        M_proc_cache_size = std::to_string( info.GetProcessorCacheSize() );
        M_proc_logical_per_physical = std::to_string( info.GetLogicalProcessorsPerPhysical() );
        M_proc_clock_frequency = std::to_string( info.GetProcessorClockFrequency() );
        M_proc_logical_cpu_number = std::to_string( info.GetNumberOfLogicalCPU() );
        M_proc_physical_cpu_number = std::to_string( info.GetNumberOfPhysicalCPU() );
        M_proc_cpu_id_support = std::to_string( info.DoesCPUSupportCPUID() );
        M_proc_apic_id = std::to_string( info.GetProcessorAPICID() );

        M_mem_virtual_total = std::to_string( info.GetTotalVirtualMemory() );
        M_mem_virtual_avail = std::to_string( info.GetAvailableVirtualMemory() );
        M_mem_physical_total = std::to_string( info.GetTotalPhysicalMemory() );
        M_mem_physical_avail = std::to_string( info.GetAvailablePhysicalMemory() );
        M_mem_host_total = std::to_string( info.GetHostMemoryTotal() );
        M_mem_host_avail = std::to_string( info.GetHostMemoryAvailable("KWSHL") );
        M_mem_proc_avail = std::to_string( info.GetProcMemoryAvailable("KWSHL", "KWSPL") );
        M_mem_host_used = std::to_string( info.GetHostMemoryUsed() );
        M_mem_proc_used = std::to_string( info.GetProcMemoryUsed() );
        M_load_avg = std::to_string( info.GetLoadAverage() );

        this->updateInfo();
    }

    ~KWSys() override = default;
};


} // Sys
} // Feel


#endif // FEELPP_KWSYS_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
