#ifndef FEELPP_HWLOC_HPP
#define FEELPP_HWLOC_HPP 1

#include <feel/feelhwsys/hwsysbase.hpp>
#include <hwloc.h>
#include <sys/utsname.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <limits.h>
#include <cstdio>
#include <cstring>
#include <string>

namespace Feel
{
namespace Sys
{

//! Wrapper class using hwloc for system information.
class HwlocSys : public HwSysBase
{
public:
    HwlocSys()
    {
        M_backend = "hwloc";

        // Initialize hwloc topology
        hwloc_topology_init(&topology_);
        hwloc_topology_load(topology_);

        // OS information via uname
        struct utsname uts;
        if (uname(&uts)==0)
        {
            M_os_name    = uts.sysname;
            M_os_release = uts.release;
            M_os_version = uts.version;
            M_os_platform= uts.machine;
        }
        
        // OS flags
    #ifdef __linux__
        M_os_is_linux   = "1";
    #else
        M_os_is_linux   = "0";
    #endif
    #ifdef __APPLE__
        M_os_is_apple   = "1";
    #else
        M_os_is_apple   = "0";
    #endif
    #ifdef _WIN32
        M_os_is_windows = "1";
    #else
        M_os_is_windows = "0";
    #endif

        // Host and domain names
        char hostbuf[HOST_NAME_MAX];
        if (gethostname(hostbuf, sizeof(hostbuf))==0)
            M_host_name = hostbuf;
        
    #ifdef __linux__
        char domainbuf[HOST_NAME_MAX];
        if (getdomainname(domainbuf, sizeof(domainbuf))==0)
            M_domain_name = domainbuf;
        else
            M_domain_name.clear();
    #else
        M_domain_name.clear();
    #endif

        // Processor topology
        int logical = hwloc_get_nbobjs_by_type(topology_, HWLOC_OBJ_PU);
        int cores   = hwloc_get_nbobjs_by_type(topology_, HWLOC_OBJ_CORE);
        int sockets = hwloc_get_nbobjs_by_type(topology_, HWLOC_OBJ_SOCKET);

        M_proc_logical_cpu_number  = std::to_string(logical);
        M_proc_physical_cpu_number = std::to_string(cores);
        if (cores>0)
            M_proc_logical_per_physical = std::to_string(logical/cores);
        else
            M_proc_logical_per_physical.clear();

        // old v1.x code — removed in hwloc v2:
        // int depth = hwloc_get_type_or_below_depth(topology_, HWLOC_OBJ_CACHE);
        // int ncaches = hwloc_get_nbobjs_by_depth(topology_, depth);
        // for (int i = 0; i < ncaches; ++i)
        // {
        //   hwloc_obj_t c = hwloc_get_obj_by_depth(topology_, depth, i);
        //   if (c->type == HWLOC_OBJ_CACHE)
        //     total_cache += c->attr->cache.size;
        // }

        // new hwloc v2-compatible version:
        size_t total_cache = 0;
        const hwloc_obj_type_t cache_levels[] = 
        {
            HWLOC_OBJ_L1CACHE,
            HWLOC_OBJ_L2CACHE,
            HWLOC_OBJ_L3CACHE,
            HWLOC_OBJ_L4CACHE
        };
        for (auto level : cache_levels)
        {
            int depth = hwloc_get_type_depth(topology_, level);
            if (depth == HWLOC_TYPE_DEPTH_UNKNOWN)
                continue;  // this level not present on the machine

            int n = hwloc_get_nbobjs_by_depth(topology_, depth);
            for (int i = 0; i < n; ++i)
            {
                hwloc_obj_t c = hwloc_get_obj_by_depth(topology_, depth, i);
                // c->attr->cache.size is the per‐object cache size in bytes
                total_cache += c->attr->cache.size;
            }
        }
        M_proc_cache_size = std::to_string(total_cache);

        // 64-bit check
        M_proc_is64bits = (sizeof(void*)>=8 ? "1" : "0");

        // Vendor string (if available)
        const char* vendor = hwloc_obj_get_info_by_name(hwloc_get_root_obj(topology_), "CPUVendor");
        M_proc_vendor_name = vendor ? vendor : std::string();

        // TODO: populate vendor_id, type_id, family_id, model_id, stepping_code, serial_number
        M_proc_vendor_id.clear();
        M_proc_type_id.clear();
        M_proc_family_id.clear();
        M_proc_model_id.clear();
        M_proc_stepping_code.clear();
        M_proc_serial_number.clear();
        M_proc_extended_name.clear();
        M_proc_cpu_id_support.clear();
        M_proc_apic_id.clear();

        // Memory and load via sysinfo and getloadavg
        struct sysinfo si;
        if (sysinfo(&si)==0)
        {
            // physical
            M_mem_physical_total = std::to_string(si.totalram);
            M_mem_physical_avail = std::to_string(si.freeram);
            // virtual = RAM + swap
            unsigned long virt_total = si.totalram + si.totalswap;
            unsigned long virt_avail = si.freeram  + si.freeswap;
            M_mem_virtual_total = std::to_string(virt_total);
            M_mem_virtual_avail = std::to_string(virt_avail);

            M_mem_host_total = M_mem_physical_total;
            M_mem_host_avail = M_mem_physical_avail;
            // proc-level usage not implemented
            M_mem_host_used.clear();
            M_mem_proc_used.clear();
            M_mem_proc_avail.clear();
        }
        
        double loads[3] = {0.0,0.0,0.0};
        if (getloadavg(loads,3)>=0)
            M_load_avg = std::to_string(loads[0]);
        else
            M_load_avg.clear();
    }

    ~HwlocSys() override
    {
        hwloc_topology_destroy(topology_);
    }

    void updateCurrentState() override
    {
        // Refresh memory and load
        struct sysinfo si;
        if (sysinfo(&si)==0)
        {
            M_mem_physical_avail = std::to_string(si.freeram);
            M_mem_host_avail     = M_mem_physical_avail;
        }
        double loads[3];
        if (getloadavg(loads,3)>=0)
            M_load_avg = std::to_string(loads[0]);
    }

private:
    hwloc_topology_t topology_;
};

} // namespace Sys
} // namespace Feel

#endif // FEELPP_HWLOC_HPP
