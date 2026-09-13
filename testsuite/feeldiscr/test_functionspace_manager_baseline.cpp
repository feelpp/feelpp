/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

   SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

   SPDX-FileCopyrightText: 2026 University of Strasbourg

   SPDX-License-Identifier: LGPL-3.0-or-later
*/

#define BOOST_TEST_MODULE test_functionspace_manager_baseline
#include <feel/feelcore/testsuite.hpp>

#include <cstdint>
#include <vector>

#include <feel/feeldiscr/functionspacebuildinstrumentation.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/unitsquare.hpp>

using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
class InstrumentationScope
{
  public:
    InstrumentationScope()
        : M_wasEnabled( FunctionSpaceBuildInstrumentation::enabled() )
    {
        FunctionSpaceBuildInstrumentation::reset();
        FunctionSpaceBuildInstrumentation::setEnabled( true );
    }

    ~InstrumentationScope()
    {
        FunctionSpaceBuildInstrumentation::reset();
        FunctionSpaceBuildInstrumentation::setEnabled( M_wasEnabled );
    }

  private:
    bool M_wasEnabled;
};
} // namespace

BOOST_AUTO_TEST_CASE( repeated_whole_mesh_factories_always_build )
{
    constexpr std::uint64_t requests = 3;
    auto mesh = unitSquare( 0.25 );
    InstrumentationScope instrumentation;

    auto checkAlwaysNew = [requests]( auto&& factory )
    {
        FunctionSpaceBuildInstrumentation::reset();

        using space_ptrtype = decltype( factory() );
        std::vector<space_ptrtype> spaces;
        spaces.reserve( requests );
        for ( std::uint64_t i = 0; i < requests; ++i )
            spaces.push_back( factory() );

        auto const counts = FunctionSpaceBuildInstrumentation::counts();
        BOOST_TEST( counts.functionSpaceConstructions == requests );
        BOOST_TEST( counts.dofTableBuilds == requests );

        for ( std::size_t i = 0; i < spaces.size(); ++i )
            for ( std::size_t j = i + 1; j < spaces.size(); ++j )
                BOOST_TEST( spaces[i].get() != spaces[j].get() );
    };

    checkAlwaysNew( [&mesh]()
                    { return Pch<2>( mesh ); } );
    checkAlwaysNew( [&mesh]()
                    { return Pdh<1>( mesh ); } );
}
