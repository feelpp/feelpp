/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "grepldeim.hpp"
#include <feel/feelmor/opusapp.hpp>
#include <feel/feelmor/crb.hpp>
#include <feel/feelmor/crbmodel.hpp>
// #include <feel/feelmor/cvgstudy.hpp>


int main( int argc, char** argv )
{
    using namespace Feel;

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options( GreplDEIM<2,2>::name() )
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(eimOptions())
                           .add(makeGreplDEIMOptions())
                           .add(podOptions())
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options( GreplDEIM<2,2>::name() )),
                           _about=makeGreplDEIMAbout( GreplDEIM<2,2>::name() ) );

    if ( ioption("cvg-study") )
    {
        // Feel::CvgStudy< GreplDEIM<2,2>, CRB, CRBModel > bench;
        // bench.run( ioption("cvg-study") );
    }
    else
    {
        Feel::OpusApp< GreplDEIM<2,2>, CRB, CRBModel > bench ;
        bench.run();
    }


}
