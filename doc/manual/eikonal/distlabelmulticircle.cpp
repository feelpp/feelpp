#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feel.hpp>

#include <feel/feells/selflabel.hpp>
#include <feel/feells/distlabels_fms.hpp>
#include <feel/feells/reinit_fms.hpp>

using namespace Feel;

int main( int argc, char** argv )
{
    boost::timer chrono;

    Feel::Environment env( argc, argv );

    const double epsilon = 1e-6;

    auto mesh = unitHypercube<FEELPP_DIM>();
    auto Xh = Pch<FEELPP_ORDER>(mesh);
    auto Xh0 = Pdh<0>(mesh);

    //ellipse parameters :
    std::vector<double> x0 = {{ 0.7, 0.15, 0.3 }};
    std::vector<double> y0 = {{ 0.6, 0.15, 0.7 }};
    std::vector<double> z0 =  {{ 0, 0, 0 }};
    if (FEELPP_DIM==3)
        z0={{ 0.7, 0.5, 0.2 }};
    std::vector<double> radius = {{ 0.1, 0.1, 0.1 }};
    std::vector<bool> hasinitlabel = {{ true, true, true }};

    const int nbcircle = x0.size();
    //const int nbcircle = 2;

    std::vector<int> label (nbcircle);
    std::iota(label.begin(), label.end(), 1);
    //std::vector<double>labelnb = {{ 1., 1., 2. }};

    // ellipse function (not exactly a distance function)
    auto phi = Xh->elementPtr();
    auto initlabel = Xh->elementPtr();
    phi->setConstant( 1e8 );

    for(int i=0;i<nbcircle; i++)
    {
        auto X0 = Px() - x0[i];
        auto Y0 = Py() - y0[i];
        auto Z0 = Pz() - z0[i];

        // multi circle function (exact distance function)

        *phi = vf::project(_space=Xh, _range=elements(mesh),
                                    _expr=min(sqrt( X0*X0 + Y0*Y0 + Z0*Z0 ) - radius[i], idv(phi) ) );

        if (hasinitlabel[i])
        {
            *initlabel = vf::project(_space=Xh, _range=elements(mesh),
                                     _expr= idv(initlabel)+label[i]*vf::chi( (X0*X0+Y0*Y0+Z0*Z0) < (radius[i]*radius[i])/2.) );
        }
    }

    // Get theoretical labels and distances
    auto L0_check = Xh->elementPtr();
    L0_check->setConstant(0.);
    auto L1_check = Xh->elementPtr();
    L1_check->setConstant(0.);
    auto L2_check = Xh->elementPtr();
    L2_check->setConstant(0.);
    for( int i = 0; i < nbcircle; ++i )
    {
        auto Xi = Px() - x0[i];
        auto Yi = Py() - y0[i];
        auto Zi = Pz() - z0[i];
        auto Ri = radius[i];

        *L0_check += vf::project(
                _space=Xh,
                _range=elements(mesh),
                _expr=(Xi*Xi + Yi*Yi + Zi*Zi < Ri*Ri) * label[i]
                );
        *L1_check += vf::project(
                _space=Xh,
                _range=elements(mesh),
                _expr=label[i]*( idv(phi) >=0 && 
                    abs((sqrt(Xi*Xi+Yi*Yi+Zi*Zi)-Ri) - idv(phi)) < epsilon )
                );

        auto distToAllExceptI = Xh->elementPtr();
        distToAllExceptI->setConstant( 1e8 );
        for( int j = 0; j < nbcircle; ++j )
        {
            if( j != i )
            {
                auto Xj = Px() - x0[j];
                auto Yj = Py() - y0[j];
                auto Zj = Pz() - z0[j];
                auto Rj = radius[j];
                *distToAllExceptI = vf::project(_space=Xh, _range=elements(mesh),
                        _expr=min(sqrt( Xj*Xj + Yj*Yj + Zj*Zj ) - Rj, idv(distToAllExceptI) ) );
            }
        }

        for( int j = 0; j < nbcircle; ++j )
        {
            if( j != i )
            {
                auto Xj = Px() - x0[j];
                auto Yj = Py() - y0[j];
                auto Zj = Pz() - z0[j];
                auto Rj = radius[j];

                *L2_check += vf::project(
                        _space=Xh,
                        _range=elements(mesh),
                        _expr=label[j] * ( 
                            abs(abs(sqrt(Xi*Xi+Yi*Yi+Zi*Zi)-Ri) - abs(idv(phi))) < epsilon &&
                            abs(abs(sqrt(Xj*Xj+Yj*Yj+Zj*Zj)-Rj) - abs(idv(distToAllExceptI))) < epsilon
                            )
                        );
            }
        }
    }

    // init label must be unchanged, selflabel will modify the label
    auto lab = Xh->elementPtr();
    *lab = *initlabel;

    auto selflabelgenerator = selfLabel(Xh, Xh0);
    selflabelgenerator->setLabel( lab );

    typedef LabelDistanceFMS<std::remove_reference<decltype(*Xh)>::type> labeldistanceFMS_type;
    auto distlabel = labeldistanceFMS_type::New( Xh );

    double timeInit = chrono.elapsed();
    Feel::cout << "initialization done in " << timeInit << " s" << std::endl;

    // interface local projection method
    chrono.restart();
    auto phiILP = vf::project(Xh, elements(mesh),
                                  idv(phi)
                                  / sqrt( inner( gradv(phi), gradv(phi) ) ) );
    double timeILP = chrono.elapsed();
    Feel::cout << "interface local projection done in " << timeILP << " s" << std::endl;

    chrono.restart();
    selflabelgenerator->updateLabel(phi);
    double timeUpdateLabel = chrono.elapsed();
    Feel::cout << "label update done in " << timeUpdateLabel << " s" << std::endl;

    chrono.restart();
    distlabel->setSelfLabel( selflabelgenerator->getLabel() );
    distlabel->run( phiILP );
    double timeDistLabel = chrono.elapsed();
    Feel::cout << "distlabel fast-marching done in " << timeDistLabel << " s" << std::endl;

    chrono.restart();
    auto fm = fms( Xh );
    auto dist = fm->march( phiILP );
    double timeFastMarching = chrono.elapsed();
    Feel::cout << "standard fast-marching done in " << timeFastMarching << " s" << std::endl;

    auto exp = exporter(_mesh=mesh, _name="distlabel");
    exp->step(0)->add("phi", *phi );
    exp->step(0)->add("initlabel", *initlabel );
    exp->step(0)->add("L0", *(selflabelgenerator->getLabel()) );
    exp->step(0)->add("L0_P0", *(selflabelgenerator->getP0Label()) );
    exp->step(0)->add("L1", *(distlabel->getNearestNeighbourLabel()) );
    exp->step(0)->add("L2", *(distlabel->getNextNearestNeighbourLabel()) );
    exp->step(0)->add("dist1", *(distlabel->getNearestNeighbourDistance()) );
    exp->step(0)->add("dist2", *(distlabel->getNextNearestNeighbourDistance()) );
    exp->step(0)->add("dist", dist );
    exp->step(0)->add("L0_check", *L0_check );
    exp->step(0)->add("L1_check", *L1_check );
    exp->step(0)->add("L2_check", *L2_check );
    exp->step(0)->addRegions();
    exp->save();

}
