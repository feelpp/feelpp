/* this file is generated automatically */
#include "stokesdeim.hpp"
#include <feel/feelcrb/opusapp.hpp>
#include <feel/feelcrb/cvgstudy.hpp>

int main( int argc, char** argv )
{
    using namespace Feel;
    po::options_description options("stokes deim options");
    options.add_options()
        ( "cvg-study", Feel::po::value<int>()->default_value( 0 ), "")
        ( "fem", Feel::po::value<bool>()->default_value( false ), "Run Fem Simulation")
        ( "mu1", Feel::po::value<double>()->default_value( 0.3 ), "mu1 for FEM")
        ( "mu2", Feel::po::value<double>()->default_value( 0.3 ), "mu2 for FEM");

    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=opusapp_options("stokesdeim")
                           .add(crbOptions())
                           .add(crbSEROptions())
                           .add(crbSaddlePointOptions())
                           .add(eimOptions())
                           .add(podOptions())
                           .add(options )
                           .add(backend_options("backend-primal"))
                           .add(backend_options("backend-dual"))
                           .add(backend_options("backend-l2"))
                           .add(bdf_options("StokesDeim")),
                           _about=makeStokesDeimAbout( "stokesdeim" ) );


    if ( ioption("cvg-study") )
    {
        Feel::CvgStudy<Feel::StokesDeim,CRBSaddlePoint,CRBModelSaddlePoint> app;
        app.run(ioption("cvg-study"));
    }
    else if ( boption("fem") )
    {
        auto model = boost::make_shared<CRBModelSaddlePoint<StokesDeim>>(crb::stage::offline);
        auto mu = model->parameterSpace()->element();
        mu[0] = doption("mu1");
        mu[1] = doption("mu2");
        auto U = model->solveFemUsingAffineDecompositionFixedPoint( mu );
        auto u = U.template element<0>();
        auto p = U.template element<1>();
        auto d = model->meshDisplacementField(mu);

        auto e = exporter( _mesh=model->functionSpace()->mesh(), _name="stokesdeimFEM" );
        e->add( "u_"+mu.toString(), u );
        e->add( "p_"+mu.toString(), p );
        e->add( "d_"+mu.toString(), *d );
        e->save();
    }
    else
    {
        Feel::OpusApp<Feel::StokesDeim,CRBSaddlePoint,CRBModelSaddlePoint> app;
        app.run();
    }

}
