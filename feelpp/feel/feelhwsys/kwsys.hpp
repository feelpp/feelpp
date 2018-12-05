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
        // Use KWSys.
        M_backend = "kwsys";

        M_sysInfo.RunCPUCheck();
        M_sysInfo.RunOSCheck();
        M_sysInfo.RunMemoryCheck();

        M_os_name = M_sysInfo.GetOSName();
        M_os_is_linux = std::to_string( M_sysInfo.GetOSIsLinux() );
        M_os_is_apple = std::to_string( M_sysInfo.GetOSIsApple() );
        M_os_is_windows = std::to_string( M_sysInfo.GetOSIsWindows() );
        M_os_release = M_sysInfo.GetOSRelease();
        M_os_version = M_sysInfo.GetOSVersion();
        M_os_platform = M_sysInfo.GetOSPlatform();
        M_host_name = M_sysInfo.GetHostname();
        M_domain_name = M_sysInfo.GetFullyQualifiedDomainName();

        M_proc_is64bits = std::to_string( M_sysInfo.Is64Bits() );
        M_proc_vendor_name = M_sysInfo.GetVendorString();
        M_proc_vendor_id = M_sysInfo.GetVendorID();
        M_proc_type_id = M_sysInfo.GetTypeID();
        M_proc_family_id = M_sysInfo.GetFamilyID();
        M_proc_model_id = M_sysInfo.GetModelID();
        M_proc_extended_name = M_sysInfo.GetExtendedProcessorName();
        M_proc_stepping_code = M_sysInfo.GetSteppingCode();
        M_proc_serial_number = M_sysInfo.GetProcessorSerialNumber();

        M_proc_cache_size = std::to_string( M_sysInfo.GetProcessorCacheSize() );
        M_proc_logical_per_physical = std::to_string( M_sysInfo.GetLogicalProcessorsPerPhysical() );
        M_proc_clock_frequency = std::to_string( M_sysInfo.GetProcessorClockFrequency() );
        M_proc_logical_cpu_number = std::to_string( M_sysInfo.GetNumberOfLogicalCPU() );
        M_proc_physical_cpu_number = std::to_string( M_sysInfo.GetNumberOfPhysicalCPU() );
        M_proc_cpu_id_support = std::to_string( M_sysInfo.DoesCPUSupportCPUID() );
        M_proc_apic_id = std::to_string( M_sysInfo.GetProcessorAPICID() );

        M_mem_virtual_total = std::to_string( M_sysInfo.GetTotalVirtualMemory() );
        M_mem_virtual_avail = std::to_string( M_sysInfo.GetAvailableVirtualMemory() );
        M_mem_physical_total = std::to_string( M_sysInfo.GetTotalPhysicalMemory() );
        M_mem_physical_avail = std::to_string( M_sysInfo.GetAvailablePhysicalMemory() );
        M_mem_host_total = std::to_string( M_sysInfo.GetHostMemoryTotal() );
        M_mem_host_avail = std::to_string( M_sysInfo.GetHostMemoryAvailable("KWSHL") );
        M_mem_proc_avail = std::to_string( M_sysInfo.GetProcMemoryAvailable("KWSHL", "KWSPL") );
        M_mem_host_used = std::to_string( M_sysInfo.GetHostMemoryUsed() );
        M_mem_proc_used = std::to_string( M_sysInfo.GetProcMemoryUsed() );
        M_load_avg = std::to_string( M_sysInfo.GetLoadAverage() );
    }

    ~KWSys() override = default;

    void updateCurrentState() override
        {
            M_mem_host_used = std::to_string( M_sysInfo.GetHostMemoryUsed() );
            M_mem_proc_used = std::to_string( M_sysInfo.GetProcMemoryUsed() );
        }
private :
    feelpp_kwsys::SystemInformation M_sysInfo;
};


} // Sys
} // Feel


#endif // FEELPP_KWSYS_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
