/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#define BOOST_TEST_MODULE modelmesh basis testsuite
#include <feel/feelcore/testsuite.hpp>

#include <string_view>

#include <feel/feelmodels/modelcore/modelmeshes.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

FEELPP_ENVIRONMENT_NO_OPTIONS

namespace
{
template <typename MeshType>
bool
hasSupportedBasis( std::string_view basis )
{
    bool found = false;
    auto const supportedBases = ModelMesh<uint32_type>::template basisFieldTypeSupported<MeshType>();
    hana::for_each( supportedBases, [&found,basis]( auto const& b )
                    {
                        if ( basis == std::string_view{ hana::at_c<0>( b ) } )
                            found = true;
                    } );
    return found;
}
}

BOOST_AUTO_TEST_SUITE( modelmesh_basis )

BOOST_AUTO_TEST_CASE( nedelec_basis_is_registered_only_for_simplex_meshes )
{
    using simplex_mesh_type = Mesh<Simplex<2>>;
    using hypercube_mesh_type = Mesh<Hypercube<2>>;

    BOOST_CHECK( hasSupportedBasis<simplex_mesh_type>( "Ned1h0" ) );
    BOOST_CHECK( !hasSupportedBasis<hypercube_mesh_type>( "Ned1h0" ) );
    BOOST_CHECK( hasSupportedBasis<hypercube_mesh_type>( "Pchv1" ) );
}

BOOST_AUTO_TEST_SUITE_END()
