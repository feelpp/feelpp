/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 */
#define BOOST_TEST_MODULE operatorinterpolation_basisexpr testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/ned1h.hpp>
#include <feel/feeldiscr/traits.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>


// FEELPP_ENVIRONMENT_WITH_OPTIONS( test_operatorinterpolation::makeAbout(),
//                                  test_operatorinterpolation::makeOptions() )
FEELPP_ENVIRONMENT_NO_OPTIONS

using namespace Feel;


template <typename Space1Type,typename Space2Type,typename RangeType>
auto executeTestBaseId( std::string const& caseName,
                        std::shared_ptr<Space1Type> Vh1, std::string const& tagDomain,
                        std::shared_ptr<Space2Type> Vh2, std::string const& tagImage,
                        RangeType const& rangeInterp )
{
    auto u1 = Vh1->element();
    auto u2 = Vh2->element();

    //auto exprDomain = -Px()*Py();
    //auto exprDomain = vec(-Py(),-Px());

    auto exprDomain = hana::eval_if( hana::bool_c<Space1Type::is_scalar>,
                                       [] { return Px()*Px(); },
                                       [] { return vec(-Py(),-Px()); } );
    //auto const& coeff_f_expr = expr( coeff_f.expr(), se );



    u1.on(_range=elements(support(Vh1)),_expr=exprDomain);

    auto opI = opInterpolation( _domainSpace=Vh1,_imageSpace=Vh2,
                                _range= rangeInterp,
                                _type= makeExprInterpolation( id(u1), nonconforming_t() ) );
    opI->matPtr()->multVector( u1, u2 );



    //auto grad_exact = vec(-Py(),-Px());
    //auto grad_exact = vec( exprDomain.template diff<1>("x"), exprDomain.template diff<1>("y") );

    std::string executeTestName = fmt::format("executeTestId_{}_{}_{}",caseName,tagDomain,tagImage);

    double errorL2 = normL2(_range=rangeInterp,_expr=idv(u2)-exprDomain);
    BOOST_MESSAGE(fmt::format("{} : errorL2={}",executeTestName,errorL2));
    BOOST_CHECK_SMALL( errorL2, 1e-10 );

    auto u2_exact = Vh2->element();
    u2_exact.on(_range=rangeInterp,_expr=exprDomain);

    double errorL2_u2exact = normL2(_range=rangeInterp,_expr=idv(u2)-idv(u2_exact) );
    BOOST_MESSAGE(fmt::format("{} : errorL2_u2exact={}",executeTestName,errorL2_u2exact));
    BOOST_CHECK_SMALL( errorL2_u2exact, 1e-10 );

    return std::make_tuple(std::move(u1),std::move(u2),std::move(u2_exact));
}




template <typename Space1Type,typename Space2Type,typename RangeType>
auto executeTestBaseGrad( std::string const& caseName,
                          std::shared_ptr<Space1Type> Vh1, std::string const& tagDomain,
                          std::shared_ptr<Space2Type> Vh2, std::string const& tagImage,
                          RangeType const& rangeInterp )
{
    auto u1 = Vh1->element();
    auto u2 = Vh2->element();

    auto exprDomain = -Px()*Py();
    //U1.setConstant(1.);
    u1.on(_range=elements(support(Vh1)),_expr=-Px()*Py());
    //u1.on(_range=elements(support(Vh1)),_expr=-Px()*Py() + 0.5*Px() - Py() );
    //u2.on(_range=elements(support(Vh2)),_expr=Px()*Py()*vec(-cst(0.5),cst(1.)));
    //u2.on(_range=elements(support(Vh2)),_expr=(1+Px()*Py())*vec(-(Py()-cst(0.5)),Px()-cst(1.)));
    //u2.setConstant(2.);

    auto opI = opInterpolation( _domainSpace=Vh1,_imageSpace=Vh2,
                                _range= rangeInterp,
                                _type= makeExprInterpolation( trans(grad(u1)), nonconforming_t() ) );
    opI->matPtr()->multVector( u1, u2 );



    auto grad_exact = vec(-Py(),-Px());
    //auto grad_exact = vec( exprDomain.template diff<1>("x"), exprDomain.template diff<1>("y") );

    std::string executeTestName = fmt::format("executeTestGrad_{}_{}_{}",caseName,tagDomain,tagImage);

    double errorL2 = normL2(_range=rangeInterp,_expr=idv(u2)-grad_exact);
    BOOST_MESSAGE(fmt::format("{} : errorL2={}",executeTestName,errorL2));
    if ( is_hcurl_conforming_v<typename Space2Type::fe_type> )
    {
        BOOST_CHECK_SMALL( errorL2, 1. );
    }
    else
    {
        BOOST_CHECK_SMALL( errorL2, 1e-10 );
    }


    auto u2_exact = Vh2->element();
    u2_exact.on(_range=rangeInterp,_expr=grad_exact);

    double errorL2_u2exact = normL2(_range=rangeInterp,_expr=idv(u2)-idv(u2_exact) );
    BOOST_MESSAGE(fmt::format("{} : errorL2_u2exact={}",executeTestName,errorL2_u2exact));
    BOOST_CHECK_SMALL( errorL2_u2exact, 1e-10 );

    return std::make_tuple(std::move(u1),std::move(u2),std::move(u2_exact));
}


BOOST_AUTO_TEST_SUITE( operatorinterpolation_basisexpr_suite )

BOOST_AUTO_TEST_CASE( test_id )
{
    using mesh_t = Mesh<Simplex<2>>;
    auto mesh = loadMesh( _mesh = new mesh_t,_filename="$cfgdir/test_operatorinterpolation_basisexpr_2d.geo" );

    auto rangeElements = elements(mesh);
    auto rangeFacesInterface = markedfaces(mesh,"Interface");
    auto rangeEltOmega1 = markedelements(mesh,"Omega1");
    auto rangeEltOmega2 = markedelements(mesh,"Omega2");
    auto submeshInterface = createSubmesh(_mesh=mesh,_range=rangeFacesInterface );
    auto submeshOmega1 = createSubmesh(_mesh=mesh,_range=rangeEltOmega1);
    auto submeshOmega2 = createSubmesh(_mesh=mesh,_range=rangeEltOmega2);

    auto rangeSubmeshInterfaceElt = elements(submeshInterface);
    auto rangeSubmeshOmega2Elt = elements(submeshOmega2);
    auto rangeSubmeshOmega2FacesInterface = markedfaces(submeshOmega2,"Interface");


    auto e = exporter( _mesh=mesh, _name="test_id_Mesh" );
    e->addRegions();

    auto eOmega2 = exporter( _mesh=submeshOmega2,_name="test_id_Omega2" );
    eOmega2->addRegions();

    auto eInterface = exporter( _mesh=submeshInterface,_name="test_id_Interface" );
    eInterface->addRegions();

    auto spacePch1 = Pch<1>(mesh);
    auto spacePch2 = Pch<2>(mesh);
    auto spacePchv1 = Pchv<1>(mesh);
    auto spacePchv2 = Pchv<2>(mesh);
    auto spaceNed1h = Ned1h<0>(mesh);

    auto spacePchv1_omega2_submesh = Pchv<1>(submeshOmega2);
    auto spaceNed1h_omega2_submesh = Ned1h<0>(submeshOmega2);

    auto spacePch1_interface_submesh = Pch<1>(submeshInterface);
    auto spacePch2_interface_submesh = Pch<2>(submeshInterface);
    auto spacePchv1_interface_submesh = Pchv<1>(submeshInterface);


    auto spacePch1_omega1_range = Pch<1>(mesh,rangeEltOmega1);
    auto spacePch2_omega1_range = Pch<2>(mesh,rangeEltOmega1);
    auto spacePchv1_omega1_range = Pchv<1>(mesh,rangeEltOmega1);
    auto spacePchv2_omega1_range = Pchv<2>(mesh,rangeEltOmega1);
    auto spacePchv1_omega2_range = Pchv<1>(mesh,rangeEltOmega2);
    auto spaceNed1h_omega2_range = Ned1h<0>(mesh,rangeEltOmega2);

    auto executeTestId = []( std::string const& caseName,
                             auto const& spaceDomain, std::string const& tagDomain,
                             auto const& spaceImage,  std::string const& tagImage,
                             auto const& rangeInterp,
                             auto const& exporterDomain, auto const& exporterImage) {
                             auto [ud,ui,ue] = executeTestBaseId(caseName,spaceDomain,tagDomain,spaceImage,tagImage,rangeInterp);
                             exporterDomain->add( fmt::format("ud_{}_{}",caseName,tagDomain), ud );
                             exporterImage->add( fmt::format("ui_{}_{}",caseName,tagImage), ui );
                             exporterImage->add( fmt::format("ue_{}_{}",caseName,tagImage), ue );
                         };


    //--------------------------------------
    // case 0 : same support
    //executeTestId( "case0_elt", spacePch1, "Pch1", spacePch1, "Pch1", rangeElements, e, e );
    //executeTest1( "case0_elt", spacePch1, "Pch1", spaceNed1h, "Ned1h", rangeElements, e, e );
    executeTestId( "case0_0_elt", spacePchv1, "Pchv1", spacePchv1, "Pchv1", rangeElements, e, e );
    executeTestId( "case0_1_elt", spacePchv1, "Pchv1", spacePchv2, "Pchv2", rangeElements, e, e );
    executeTestId( "case0_2_elt", spacePchv1, "Pchv2", spacePchv2, "Pchv1", rangeElements, e, e );
    executeTestId( "case0_3_elt", spacePchv1, "Pchv2", spacePchv2, "Pchv2", rangeElements, e, e );

    executeTestId( "case0_0_face_interface", spacePchv1, "Pchv1", spacePchv1, "Pchv1", rangeFacesInterface, e, e );
    executeTestId( "case0_1_face_interface", spacePchv1, "Pchv1", spacePchv2, "Pchv2", rangeFacesInterface, e, e );
    executeTestId( "case0_2_face_interface", spacePchv1, "Pchv2", spacePchv2, "Pchv1", rangeFacesInterface, e, e );
    executeTestId( "case0_3_face_interface", spacePchv1, "Pchv2", spacePchv2, "Pchv2", rangeFacesInterface, e, e );

    //--------------------------------------
    // case 1 : supports intersection = Omega2
    // use submesh
    executeTestId( "case1_0_elt_sm", spacePchv1, "Pchv1", spacePchv1_omega2_submesh, "Pchv1", rangeSubmeshOmega2Elt/*rangeEltOmega2*/, e, eOmega2 );

    //--------------------------------------
    // case 2 : image submesh = interface
    executeTestId( "case2_elt_ism", spacePchv1, "Pchv1", spacePchv1_interface_submesh, "Pchv1", rangeSubmeshInterfaceElt, e, eInterface );
    if ( Environment::numberOfProcessors() == 1  ) // TO INVESTIGATE
        executeTestId( "case2_face_interface_dsm", spacePchv1_interface_submesh, "Pchv1", spacePchv1, "Pchv1", rangeFacesInterface, eInterface, e );
    if ( Environment::numberOfProcessors() == 1  )
        executeTestId( "case2_face_interface_range_dsm", spacePchv1_interface_submesh, "Pchv1", spacePchv1_omega2_range, "Pchv1", rangeFacesInterface, eInterface, e );
    executeTestId( "case2_elt_dsm_ism", spacePchv1_interface_submesh, "Pchv1", spacePchv1_interface_submesh, "Pchv1", rangeSubmeshInterfaceElt, eInterface, eInterface );

    //--------------------------------------
    // case 3 : supports disjoint , intersection = interface
    if ( Environment::numberOfProcessors() == 1  )
        executeTestId( "case3_face_interface", spacePchv1_omega1_range, "Pchv1", spacePchv1_omega2_range, "Pchv1", rangeFacesInterface, e, e );


    e->save();
    eOmega2->save();
    eInterface->save();


}
#if 1
BOOST_AUTO_TEST_CASE( test_grad )
{
    using mesh_t = Mesh<Simplex<2>>;
    auto mesh = loadMesh( _mesh = new mesh_t,_filename="$cfgdir/test_operatorinterpolation_basisexpr_2d.geo" );

    auto rangeElements = elements(mesh);
    auto rangeFacesInterface = markedfaces(mesh,"Interface");
    auto rangeEltOmega1 = markedelements(mesh,"Omega1");
    auto rangeEltOmega2 = markedelements(mesh,"Omega2");
    auto submeshInterface = createSubmesh(_mesh=mesh,_range=rangeFacesInterface );
    auto submeshOmega1 = createSubmesh(_mesh=mesh,_range=rangeEltOmega1);
    auto submeshOmega2 = createSubmesh(_mesh=mesh,_range=rangeEltOmega2);

    auto rangeSubmeshInterfaceElt = elements(submeshInterface);
    auto rangeSubmeshOmega2Elt = elements(submeshOmega2);
    auto rangeSubmeshOmega2FacesInterface = markedfaces(submeshOmega2,"Interface");


    auto e = exporter( _mesh=mesh, _name="test_grad_Mesh" );
    e->addRegions();

    auto eOmega2 = exporter( _mesh=submeshOmega2,_name="test_grad_Omega2" );
    eOmega2->addRegions();

    auto eInterface = exporter( _mesh=submeshInterface,_name="test_grad_Interface" );
    eInterface->addRegions();

    auto spacePch1 = Pch<1>(mesh);
    auto spacePch2 = Pch<2>(mesh);
    auto spacePchv1 = Pchv<1>(mesh);
    auto spaceNed1h = Ned1h<0>(mesh);

    auto spacePchv1_omega2_submesh = Pchv<1>(submeshOmega2);
    auto spaceNed1h_omega2_submesh = Ned1h<0>(submeshOmega2);

    auto spacePch1_interface_submesh = Pch<1>(submeshInterface);
    auto spacePch2_interface_submesh = Pch<2>(submeshInterface);
    auto spacePchv1_interface_submesh = Pchv<1>(submeshInterface);


    auto spacePch1_omega1_range = Pch<1>(mesh,rangeEltOmega1);
    auto spacePch2_omega1_range = Pch<2>(mesh,rangeEltOmega1);
    auto spacePchv1_omega2_range = Pchv<1>(mesh,rangeEltOmega2);
    auto spaceNed1h_omega2_range = Ned1h<0>(mesh,rangeEltOmega2);


    auto executeTestGrad = []( std::string const& caseName,
                               auto const& spaceDomain, std::string const& tagDomain,
                               auto const& spaceImage,  std::string const& tagImage,
                               auto const& rangeInterp,
                               auto const& exporterDomain, auto const& exporterImage) {
                               auto [ud,ui,ue] = executeTestBaseGrad(caseName,spaceDomain,tagDomain,spaceImage,tagImage,rangeInterp);
                               exporterDomain->add( fmt::format("ud_{}_{}",caseName,tagDomain), ud );
                               exporterImage->add( fmt::format("ui_{}_{}",caseName,tagImage), ui );
                               exporterImage->add( fmt::format("ue_{}_{}",caseName,tagImage), ue );
                           };


    //--------------------------------------
    // case 0 : same support
    executeTestGrad( "case0_elt", spacePch1, "Pch1", spaceNed1h, "Ned1h", rangeElements, e, e );
    executeTestGrad( "case0_elt", spacePch2, "Pch2", spacePchv1, "Pchv1", rangeElements, e, e );
    executeTestGrad( "case0_face_interface", spacePch1, "Pch1", spaceNed1h, "Ned1h", rangeFacesInterface, e, e );
    executeTestGrad( "case0_face_interface", spacePch2, "Pch2", spacePchv1, "Pchv1", rangeFacesInterface, e, e );

    //--------------------------------------
    // case 1 : supports intersection = Omega2
    // use submesh
    executeTestGrad( "case1_elt_sm", spacePch1, "Pch1", spaceNed1h_omega2_submesh, "Ned1h", rangeSubmeshOmega2Elt/*rangeEltOmega2*/, e, eOmega2 );
    executeTestGrad( "case1_elt_sm", spacePch2, "Pch2", spacePchv1_omega2_submesh, "Pchv1", rangeSubmeshOmega2Elt/*rangeEltOmega2*/, e, eOmega2 );
    if ( Environment::numberOfProcessors() == 1  )
        executeTestGrad( "case1_face_interface_sm", spacePch1, "Pch1", spaceNed1h_omega2_submesh, "Ned1h", rangeSubmeshOmega2FacesInterface, e, eOmega2 );
    executeTestGrad( "case1_face_interface_sm", spacePch2, "Pch2", spacePchv1_omega2_submesh, "Pchv1", rangeSubmeshOmega2FacesInterface, e, eOmega2 );

    // use space range
    if ( Environment::numberOfProcessors() == 1  )
        executeTestGrad( "case1_elt_range", spacePch1, "Pch1", spaceNed1h_omega2_range, "Ned1h", rangeEltOmega2, e, e );
    executeTestGrad( "case1_elt_range", spacePch2, "Pch2", spacePchv1_omega2_range, "Pchv1", rangeEltOmega2, e, e );
    if ( Environment::numberOfProcessors() == 1  )
        executeTestGrad( "case1_face_interface_range", spacePch1, "Pch1", spaceNed1h_omega2_range, "Ned1h", rangeFacesInterface, e, e );
    if ( Environment::numberOfProcessors() == 1  )
        executeTestGrad( "case1_face_interface_range", spacePch2, "Pch2", spacePchv1_omega2_range, "Pchv1", rangeFacesInterface, e, e );

    //--------------------------------------
    // case 2 : image submesh = interface
    //executeTestGrad( "case1_elt", spacePch1, "Pch1", spacePchv1_interface_submesh, "Ned1h", rangeSubmeshInterfaceElt, e, eInterface );
    executeTestGrad( "case2_elt_ism", spacePch2, "Pch2", spacePchv1_interface_submesh, "Pchv1", rangeSubmeshInterfaceElt, e, eInterface );
#if 0
    // TODO FIX THESE TESTS
    executeTestGrad( "case2_face_interface_dsm", spacePch2_interface_submesh, "Pch2", spacePchv1, "Pchv1", rangeFacesInterface, eInterface, e );
    executeTestGrad( "case2_face_interface_range_dsm", spacePch2_interface_submesh, "Pch2", spacePchv1_omega2_range, "Pchv1", rangeFacesInterface, eInterface, e );
    //executeTestGrad( "case2_elt_dsm_ism", spacePch2_interface_submesh, "Pch2", spacePchv1_interface_submesh, "Pchv1", rangeSubmeshInterfaceElt, eInterface, eInterface );
#endif

    //--------------------------------------
    // case 3 : supports disjoint , intersection = interface
    if ( Environment::numberOfProcessors() == 1  )
    {
        executeTestGrad( "case3_face_interface", spacePch1_omega1_range, "Pch1", spaceNed1h_omega2_range, "Ned1h", rangeFacesInterface, e, e );
        executeTestGrad( "case3_face_interface", spacePch2_omega1_range, "Pch2", spacePchv1_omega2_range, "Pchv1", rangeFacesInterface, e, e );
    }


    e->save();
    eOmega2->save();
    eInterface->save();
}
#endif

BOOST_AUTO_TEST_SUITE_END()
