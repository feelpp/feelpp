// tests/test_hwlocsys_boostut.cpp

#include <boost/ut.hpp>
#include <feel/feelhwsys/kwsys.hpp>
#include <string>

using namespace boost::ut;

int main()
{
    "HwlocSys basic self-test"_test = []() {
        Feel::Sys::HwlocSys sys;

        // 1) At least one logical CPU
        int logical = std::stoi(sys.procLogicalCpuNumber());
        expect(logical > 0_i)
            << "Expected at least 1 logical CPU, got " << logical;

        // 2) At least one physical core
        int physical = std::stoi(sys.procPhysicalCpuNumber());
        expect(physical > 0_i)
            << "Expected at least 1 physical core, got " << physical;

        // 3) Cache size should be non-zero
        unsigned long long cache = std::stoull(sys.procCacheSize());
        expect(cache > 0_ull)
            << "Expected total cache > 0 bytes, got " << cache;

        // 4) Total RAM should be non-zero
        unsigned long long ram = std::stoull(sys.memPhysicalTotal());
        expect(ram > 0_ull)
            << "Expected total physical RAM > 0, got " << ram;

        // 5) Load average should be ≥ 0.0
        double load = std::stod(sys.loadAvg());
        expect(load >= 0.0_d)
            << "Expected non-negative load average, got " << load;
    };
}