/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-02-20

  Copyright (C) 2013 Université de Strasbourg

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_env.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-02-20
 */
#define USE_BOOST_TEST 1
// Boost.Test

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE kokkos testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/kokkos.hpp>
#include <hwloc.h>
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description integrationoptions( "Test Kokkos options" );
    integrationoptions.add_options()
        ( "N", Feel::po::value<int>()->default_value( 20 ), "N value" )
    ;
    return integrationoptions.add( Feel::feel_options() );
}

inline Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_kokkos",
                           "test_kokkos",
                           "0.2",
                           "Kokkos tests",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2024 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

template <typename ExecSpace>
void test_kokkos(const std::string& label, int N) 
{
    // Define Views in the execution space's memory space
    using MemorySpace = typename ExecSpace::memory_space;
    Kokkos::View<double*, MemorySpace> x("x", N);
    Kokkos::View<double*, MemorySpace> y("y", N);
    Kokkos::View<double*, MemorySpace> z("z", N);

    // Initialize x and y
    Kokkos::parallel_for("Initialize", Kokkos::RangePolicy<ExecSpace>(0, N), KOKKOS_LAMBDA(int i) {
        x(i) = 1.0;
        y(i) = 2.0;
    });

    // Fence to ensure initialization is complete
    Kokkos::fence();

    // Start timing
    Kokkos::Timer timer;

    // Perform z = x + y
    Kokkos::parallel_for(label, Kokkos::RangePolicy<ExecSpace>(0, N), KOKKOS_LAMBDA(int i) {
        z(i) = x(i) + y(i);
    });

    // Fence to ensure computation is complete
    Kokkos::fence();

    // Stop timing
    double elapsed_time = timer.seconds();

    // Output the timing result
    std::cout << fmt::format( "{} | {} execution time : {} seconds", N, label, elapsed_time ) << std::endl;

    bool success = true;
    Kokkos::parallel_reduce( "Verify", Kokkos::RangePolicy<ExecSpace>(0, N),
                             KOKKOS_LAMBDA(const int i, bool& local_success) {
                                if (z(i)!=3.0) 
                                {
                                    local_success = false;
                                }
                             },
                            Kokkos::LAnd<bool>(success)  // Logical AND reduction
                            );
    BOOST_CHECK_MESSAGE(success, fmt::format("{} computation failed",label));
}

template <typename ExecSpace, typename HostMemSpace = Kokkos::HostSpace>
void test_kokkos_2(const std::string& label, int N)
{
    using MemorySpace = typename ExecSpace::memory_space;
    using namespace Feel;
    using ValueType = double;
    constexpr size_t dataSize = 1024 * 1024; // 1 million elements

    int rank = Environment::rank();
    int size = Environment::numberOfProcessors();

    int num_gpus = 3;
    //hipGetDeviceCount(&num_gpus);
    int device_id = rank % num_gpus;

    // Create a Kokkos View on the device
    Kokkos::View<ValueType*, MemorySpace> d_data("d_data", dataSize);

    // Perform computation on the GPU: Initialize data
    Kokkos::parallel_for("InitData", Kokkos::RangePolicy<ExecSpace>(0, dataSize),
        KOKKOS_LAMBDA(const size_t i) {
            d_data(i) = static_cast<ValueType>(i);
    });

    // Create a host mirror of the device data
    auto h_data = Kokkos::create_mirror_view(d_data);

    // Copy data from device to host
    Kokkos::deep_copy(h_data, d_data);

    // Now use Boost.MPI to communicate the data
    if (rank == 0) 
    {
        // Process 0 sends data to all other processes
        for (int dest = 1; dest < size; ++dest) 
        {
            // Send data using Boost.MPI
            Environment::worldComm().send(dest, 0, h_data.data(), dataSize);
        }
    } 
    else 
    {
        // Other processes receive data from process 0
        Environment::worldComm().recv(0, 0, h_data.data(), dataSize);

        // Copy data back to device
        Kokkos::deep_copy(d_data, h_data);

        // Perform further computations on the GPU
        Kokkos::parallel_for("ProcessData", Kokkos::RangePolicy<ExecSpace>(0, dataSize),
                            KOKKOS_LAMBDA(const size_t i) 
                            {
                                d_data(i) *= 2.0;
                            });

        // Copy results back to host for verification
        Kokkos::deep_copy(h_data, d_data);

        // Output some data for verification
        if (rank == 1) 
        {
            std::cout << "Process " << rank << " first element: " << h_data(0) << std::endl;
        }
    }
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )


BOOST_AUTO_TEST_SUITE( kokkos )

BOOST_AUTO_TEST_CASE( kokkos_1 )
{
    using namespace Feel;
    const int N = std::pow<int>(2,ioption(_name="N"));
    LOG(INFO) << "N = " << N << std::endl;

    // Run computation in Serial execution space
    test_kokkos<Kokkos::Serial>("Serial", N);

    // Run computation in Threads execution space
    test_kokkos<Kokkos::Threads>("Threads", N);

    // Run computation in HIP execution space (if enabled)
#ifdef KOKKOS_ENABLE_HIP
    test_kokkos<Kokkos::HIP>("HIP", N);
#endif
}

BOOST_AUTO_TEST_CASE( kokkos_2 )
{
    using namespace Feel;
    const int N = std::pow<int>(2,ioption(_name="N"));
    LOG(INFO) << "N = " << N << std::endl;

#ifdef KOKKOS_ENABLE_HIP
    test_kokkos_2<Kokkos::HIP,Kokkos::HIPSpace>("HIP", N);
#endif

}

BOOST_AUTO_TEST_SUITE_END()


