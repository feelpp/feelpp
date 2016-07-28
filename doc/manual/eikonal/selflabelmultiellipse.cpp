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
    Feel::Environment env( argc, argv );

    auto mesh = unitHypercube<FEELPP_DIM>();
    auto Xh = Pch<FEELPP_ORDER>(mesh);
    auto Xh0 = Pdh<0>(mesh);

    //ellipse parameters :
    std::vector<double>x0 = {{ 0.7, 0.15, 0,3 }};
    std::vector<double>y0 = {{ 0.6, 0.15, 0.7 }};
    std::vector<double>z0 =  {{0,0,0}};
    if (FEELPP_DIM==3)
        z0={{ 0.7, 0.5, 0.2 }};
    std::vector<double>aAxis = {{ 0.1, 0.1, 0.1 }};
    std::vector<double>bAxis = {{ 0.3, 0.1, 0.15 }};
    std::vector<bool>hasinitlabel = {{ true, true, false }};
    //std::vector<double>labelnb = {{ 1., 1., 2. }};

    auto ellipseShape = Xh->elementPtr();
    auto initlabel = Xh->elementPtr();
    ellipseShape->setConstant( 1e8 );
    const int nbellipse = x0.size();

    for(int i=0;i<nbellipse; i++)
    {
        auto X0 = Px() - x0[i];
        auto Y0 = Py() - y0[i];
        auto Z0 = Pz() - z0[i];

        // ellipse function (not exactly a distance function)

        *ellipseShape = vf::project(_space=Xh, _range=elements(mesh),
                                    _expr=min(sqrt( (X0/aAxis[i]) * (X0/aAxis[i])
                                                    + (Y0/bAxis[i]) * (Y0/bAxis[i])
                                                    + (Z0/bAxis[i]) * (Z0/bAxis[i]) ) - 1, idv(ellipseShape) ) );

        if (hasinitlabel[i])
        {
            *initlabel = vf::project(_space=Xh, _range=elements(mesh),
                                     _expr= idv(initlabel)+(i+1)*vf::chi( (X0*X0+Y0*Y0+Z0*Z0) < (aAxis[i]*aAxis[i])/2.) );
        }
    }

    // init label must be unchanged, selflabel will modify the label
    auto lab = Xh->elementPtr();
    *lab = *initlabel;

    auto selflabelgenerator = selfLabel(Xh, Xh0);
    selflabelgenerator->setLabel( lab );
    selflabelgenerator->updateLabel(ellipseShape);

    auto exp = exporter(_mesh=mesh, _name="selflabel");
    exp->step(0)->add("L0", *(selflabelgenerator->getLabel()) );
    exp->step(0)->add("ellipseShape", *ellipseShape );
    exp->step(0)->add("L0_P0", *(selflabelgenerator->getP0Label()) );
    exp->step(0)->add("initlabel", *initlabel );
    exp->save();
}
