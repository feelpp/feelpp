/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "grepldeim.hpp"
#include <feel/feelmor/opusapp.hpp>
#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/crbmodel.hpp>
#include <feel/feelmor/cvgstudy.hpp>

inline po::options_description makeOptions()
{
    po::options_description grepldeimoptions( "GREPL DEIM options");
    grepldeimoptions.add_options()
        ( "trainset-deim-size", Feel::po::value<int>()->default_value( 40 ), "EIM trainset is built using a equidistributed grid 40 * 40 by default")
        ( "gamma", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary Dirichlet formulation" )
        ( "cvg-study", Feel::po::value<int>()->default_value( 0 ), "size of the sampling for the convergence study. If 0 : no cvg study" )
        ;
   return grepldeimoptions.add( feel_options());
}

inline AboutData
makeAbout( std::string const& app_name)
{
    AboutData about( app_name.c_str() );
    return about;
}

int main( int argc, char** argv )
{
    using namespace Feel;

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options( GreplDEIM<2,2>::name() )
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(eimOptions())
                           .add(makeOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options( GreplDEIM<2,2>::name() )),
                           _about=makeAbout( GreplDEIM<2,2>::name() ) );

    if ( ioption("cvg-study") )
    {
        Feel::CvgStudy< GreplDEIM<2,2>, CRB, CRBModel > bench;
        bench.run( ioption("cvg-study") );
    }
    else
    {
        Feel::OpusApp< GreplDEIM<2,2>, CRB, CRBModel > bench ;
        bench.run();
    }


}
