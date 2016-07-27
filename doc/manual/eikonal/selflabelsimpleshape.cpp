#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feel.hpp>

#include <feel/feells/selflabel.hpp>

using namespace Feel;

int main( int argc, char** argv )
{
    const int dim = 2;

    Feel::Environment env( argc, argv );

    auto mesh = unitHypercube<dim>();
    auto Xh = Pch<1>(mesh);
    auto Xh0 = Pdh<0>(mesh);

        // some ellipse parameters:
    const double x0=0.6;
    const double y0=0.5;
    const double z0 = dim==2 ? 0 : 0.5;
    const double aAxis = 0.1;
    const double bAxis = 0.3;

    auto X0 = Px() - x0;
    auto Y0 = Py() - y0;
    auto Z0 = Pz() - z0;

    // ellipse function (not exactly a distance function)
    auto ellipseShape = Xh->elementPtr();
    *ellipseShape = vf::project(_space=Xh, _range=elements(mesh),
                                    _expr=sqrt( (X0/aAxis) * (X0/aAxis)
                                                + (Y0/bAxis) * (Y0/bAxis)
                                                + (Z0/bAxis) * (Z0/bAxis) ) - 1 );

    auto initlabel = Xh->elementPtr();
    *initlabel = vf::project(_space=Xh, _range=elements(mesh),
                                 _expr= vf::chi( (X0*X0+Y0*Y0+Z0*Z0)< (aAxis*aAxis)/8.) );

    auto selflabelgenerator = selfLabel(Xh, Xh0);
    selflabelgenerator->setLabel( initlabel );
    selflabelgenerator->updateLabel(ellipseShape);

    auto exp = exporter(_mesh=mesh, _name="selflabel");
    exp->step(0)->add("L0", *(selflabelgenerator->getLabel()) );
    exp->step(0)->add("L0_P0", *(selflabelgenerator->getP0Label()) );
    exp->step(0)->add("initlabel", *initlabel );
    exp->save();

}
