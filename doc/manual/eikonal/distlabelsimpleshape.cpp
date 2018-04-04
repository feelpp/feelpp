#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feel.hpp>

#include <feel/feells/selflabel.hpp>
#include <feel/feells/distlabels_fms.hpp>

using namespace Feel;

int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    auto mesh = unitHypercube<FEELPP_DIM>();
    auto Xh = Pch<FEELPP_ORDER>(mesh);
    auto Xh0 = Pdh<0>(mesh);

        // some ellipse parameters:
    const double x0=0.6;
    const double y0=0.5;
    const double z0 = FEELPP_DIM==2 ? 0 : 0.5;
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
                                                + (Z0/bAxis) * (Z0/bAxis) ) - 1. );

    auto initlabel = Xh->elementPtr();
    *initlabel = vf::project(_space=Xh, _range=elements(mesh),
                                 _expr= vf::chi( (X0*X0+Y0*Y0+Z0*Z0) < (aAxis*aAxis)/2.) );

    // init label must be unchanged, selflabel will modify the label
    auto lab = Xh->elementPtr();
    *lab = *initlabel;

    auto selflabelgenerator = selfLabel(Xh, Xh0);
    selflabelgenerator->setLabel( lab );

    typedef LabelDistanceFMS<std::remove_reference<decltype(*Xh)>::type> labeldistanceFMS_type;
    auto distlabel = labeldistanceFMS_type::New( Xh );

    selflabelgenerator->updateLabel(ellipseShape);
    distlabel->setSelfLabel( selflabelgenerator->getLabel() );
    distlabel->run( *ellipseShape );

    auto exp = exporter(_mesh=mesh, _name="distlabel");
    exp->step(0)->add("L0", *(selflabelgenerator->getLabel()) );
    exp->step(0)->add("ellipseShape", *ellipseShape );
    exp->step(0)->add("L0_P0", *(selflabelgenerator->getP0Label()) );
    exp->step(0)->add("initlabel", *initlabel );
    exp->step(0)->add("L1", *(distlabel->getNearestNeighbourLabel()) );
    exp->step(0)->add("dist1", *(distlabel->getNearestNeighbourDistance()) );
    exp->step(0)->add("L2", *(distlabel->getNextNearestNeighbourLabel()) );
    exp->step(0)->add("dist2", *(distlabel->getNextNearestNeighbourDistance()) );
    exp->step(0)->addRegions();
    exp->save();

}
