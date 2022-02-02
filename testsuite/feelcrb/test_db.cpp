/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Alexandre Ancel <alexandre.ancel@cemosis.fr>
       Date: 2016-04-06

  Copyright (C) 2016 Universite de Strasbourg

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
   \file test_db.cpp
   \author Alexandre Ancel <alexandre.ancel@cemosis.fr>
   \date 2016-04-06
 */

#define BOOST_TEST_MODULE crb db testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/unitsquare.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelcrb/opusapp.hpp>

namespace Feel
{
inline
po::options_description
makeOptions()
{
    po::options_description dboptions( "test_db options" );
    dboptions
        .add(opusapp_options("feelpp_test_db"))
        .add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"));
    /*
    dboptions.add_options()
    ( "cvg-study" , po::value<bool>()->default_value( false ), "run a convergence study if true" )
    ;
    */
    return dboptions;
}

inline
AboutData
makeAbout()
{
    AboutData about( "feelpp_test_db" ,
                     "feelpp_test_db" ,
                     "0.1",
                     "Database tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Universite de Strasbourg" );

    about.addAuthor( "Alexandre Ancel", "developer", "alexandre.ancel@cemosis.fr", "" );
    return about;

}



class TestDbHeat1d : public ModelCrbBase<ParameterSpaceX, decltype(Pch<5>(Mesh<Simplex<1>>::New()))>
{
public:

    typedef ModelCrbBase<ParameterSpaceX, decltype(Pch<5>(Mesh<Simplex<1>>::New()))> super_type;

    TestDbHeat1d()
        :
        super_type( (boost::format("TestDbHeat1d_%1%_np%2%")%soption(_name="crb.db.format") %Environment::worldComm().size() ).str() )
        {
            this->setId( boost::uuids::nil_uuid() );
        }

    //! initialisation of the model
    void initModel()
        {
            this->setFunctionSpaces( Pch<5>( loadMesh( _mesh=new Mesh<Simplex<1>> ) ) );

            Dmu->setDimension( 4 );
            //static const int N = 2;
            auto mu_min = Dmu->element();
            mu_min << 0.2, 0.2, 0.01, 0.1;
            Dmu->setMin( mu_min );
            auto mu_max = Dmu->element();
            mu_max << 50, 50, 5, 5;
            Dmu->setMax( mu_max );

            auto u = Xh->element();
            auto v = Xh->element();
            auto mesh = Xh->mesh();
            //lhs
            auto a0 = form2( _trial=Xh, _test=Xh);
            a0 = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
                integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) );
            this->addLhs( { a0 , "1" } );

            auto a1 = form2( _trial=Xh, _test=Xh);
            a1 = integrate( _range=markedelements( mesh,"k1_1" ),  _expr=gradt( u )*trans( grad( v ) )  );
            this->addLhs( { a1 , "mu0" } );

            auto a2 = form2( _trial=Xh, _test=Xh);
            a2 = integrate( _range=markedelements( mesh,"k2_1" ),  _expr=gradt( u )*trans( grad( v ) )  );
            this->addLhs( { a2 , "mu1" } );

            //rhs
            auto f0 = form1( _test=Xh );
            f0 = integrate( _range=markedfaces( mesh,"left" ), _expr=id( v ) );
            this->addRhs( { f0, "mu2" } );
            auto f1 = form1( _test=Xh );
            f1 =  integrate( _range=elements( mesh ), _expr=id( v ) );
            this->addRhs( { f1, "mu3" } );

            //output
            auto out = form1( _test=Xh );
            out = integrate( _range=markedelements( mesh,"k1_2" ), _expr=id( v )/0.2 ) +
                integrate( _range=markedelements( mesh,"k2_1" ), _expr=id( v )/0.2 );
            this->addOutput( { out, "1" } );

            auto energy = form2( _trial=Xh, _test=Xh);
            energy = integrate( _range=elements( mesh ), _expr=0.1*( gradt( u )*trans( grad( v ) ) ) ) +
                integrate( _range=markedfaces( mesh,"right" ), _expr=idt( u )*id( v ) ) +
                integrate( _range=markedelements( mesh,"k1_1" ),  _expr=0.2 * gradt( u )*trans( grad( v ) ) )  +
                integrate( _range=markedelements( mesh,"k2_1" ),  _expr=0.2 * gradt( u )*trans( grad( v ) ) )  ;
            this->addEnergyMatrix( energy );
        }

    beta_vector_light_type beta;
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false)
        {
            auto mesh = Xh->mesh();
            double output=0;
            // right hand side (compliant)
            if ( output_index == 0 )
            {
                output  = integrate( _range=markedfaces( mesh,"left" ), _expr=mu(2)*idv( u ) ).evaluate()( 0,0 );
                output += integrate( _range=elements( mesh ), _expr=mu(3)*idv( u ) ).evaluate()( 0,0 );
            }
            // output
            else if ( output_index == 1 )
            {
                output = integrate( _range=elements( mesh ),
                                    _expr=chi( ( Px() >= -0.1 ) && ( Px() <= 0.1 ) )*idv( u ) ).evaluate()( 0,0 )/0.2;
            }
            else
                throw std::logic_error( "[Heat1d::output] error with output_index : only 0 or 1 " );
            return output;
        }
};




};

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

/* This test is based on the Heat1d application */
/* It tests different setup for the databases: */
/* - Whether we create or load a new db */
/* - Whether we use the boost or hdf5 backend */
BOOST_AUTO_TEST_SUITE( crb_db )

BOOST_AUTO_TEST_CASE( crd_db_test1 )
{
    std::vector<int> wnSize;
    std::vector<int> wnduSize;

    int sampleSize = 4;
    std::vector<double> wnSample;

    // Setup initial options
    Environment::setOptionValue("crb.results-repo-name", std::string("test_db"));
    Environment::setOptionValue("crb.rebuild-database", true);

    // Execute the Heat1d app with both boost and hdf5 databases
    Environment::setOptionValue("crb.db.format", std::string("boost"));
    Feel::OpusApp<Feel::TestDbHeat1d > * app = new Feel::OpusApp<Feel::TestDbHeat1d>();
    app->run();

    auto crb = app->getCRB();
    auto WN = crb->wn();
    auto WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0]->size() > sampleSize)
    {
        auto const& wn0 = unwrap_ptr(WN[0]);
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(wn0[i]); }
    }

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::TestDbHeat1d >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0]->size() > sampleSize)
    {
        auto const& wn0 = unwrap_ptr(WN[0]);
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(wn0[i]); }
    }

    delete app;

    // Execute the Heat1d app with both boost and hdf5 databases
    // Except we are reloading the databases here
    Environment::setOptionValue("crb.rebuild-database", false);

    Environment::setOptionValue("crb.db.format", std::string("boost"));
    app = new Feel::OpusApp<Feel::TestDbHeat1d >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0]->size() > sampleSize)
    {
        auto const& wn0 = unwrap_ptr(WN[0]);
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(wn0[i]); }
    }

    delete app;

    Environment::setOptionValue("crb.db.format", std::string("hdf5"));
    app = new Feel::OpusApp<Feel::TestDbHeat1d >();
    app->run();

    crb = app->getCRB();
    WN = crb->wn();
    WNdu = crb->wndu();
    wnSize.push_back(WN.size());
    wnduSize.push_back(WNdu.size());

    if(WN.size() > 0 && WN[0]->size() > sampleSize)
    {
        auto const& wn0 = unwrap_ptr(WN[0]);
        for(int i = 0; i < sampleSize; i++)
        { wnSample.push_back(wn0[i]); }
    }

    /*
    std::cout << wnSample.size() << std::endl;

    for(int i = 0; i < wnSample.size(); i++)
    {
        if(i != 0 && (i % sampleSize) == 0)
        {
            std::cout << std::endl;
        }
        std::cout << wnSample.at(i) << " ";
    }
    std::cout << std::endl;
    */

    delete app;

    /* Check that we have the correct basis sizes */
    std::cout << "Checking basis sizes ..." << std::endl;
    for(int i = 1; i < wnSize.size(); i++)
    {
        std::cout << "WN0:size=" << wnSize[0] << " ; WN" << i << ":size=" << wnSize[i] << std::endl;
        BOOST_CHECK_MESSAGE(std::abs(wnSize[0] - wnSize[i]) <= 1, "Basis sizes OK" );
    }
    std::cout << "Basis sizes OK" << std::endl;

    /* Check that we have the correct amount of samples */
    std::cout << "Checking Number of samples ..." << std::endl;
    BOOST_CHECK(wnSample.size() == 4 * sampleSize);
    std::cout << "Number of samples OK" << std::endl;

    double tolerance = 1e-5;
    int nbGroups = wnSample.size() / sampleSize;
    std::cout << "Checking Samples ..." << std::endl;
    for(int i = 0; i < sampleSize; i++)
    {
        std::cout << "wnSample[i]:" << wnSample[i] << " ; " << wnSample[i + sampleSize] << std::endl;
        // std tolerance : std::fabs(a - b) < std::numeric_limits<double>::epsilon();
        //if(std::fabs(wnSample[i] - wnSample[i + sampleSize]) > tolerance)
        BOOST_CHECK(std::fabs(wnSample[i] - wnSample[i + sampleSize]) < tolerance);
        std::cout << "wnSample[i]:" << wnSample[i] << " ; wnSample[i + 2 * sampleSize]):" << wnSample[i + 2 * sampleSize] << std::endl;
        BOOST_CHECK(std::fabs(wnSample[i] - wnSample[i + 2 * sampleSize]) < tolerance);
        std::cout << "wnSample[i]:" << wnSample[i] << " ; wnSample[i + 3 * sampleSize]):" << wnSample[i + 3 * sampleSize] << std::endl;
        BOOST_CHECK(std::fabs(wnSample[i] - wnSample[i + 3 * sampleSize]) < tolerance);
    }
    std::cout << "Samples OK" << std::endl;

}
BOOST_AUTO_TEST_SUITE_END()
