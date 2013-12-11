#include <feel/feel.hpp>
#include <feel/feelpde/reinit_fms.hpp>

#define DIM 2


inline
Feel::AboutData
makeAboutFSI()
{
    Feel::AboutData about( "DistanceToWallsFmOnly",
                           "distanceToWallsOnly",
                           "0.1",
                           "get the distance to walls thanks to level set framework",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2013 Universite de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Vincent Doyeux", "developer", "vincent.doyeux@gmail.com", "");
    return about;
}


using namespace Feel;
using namespace Feel::vf;

void run()
{
    typedef Mesh< Simplex<DIM> > mesh_type;

    auto mesh = loadMesh( _mesh=new mesh_type );

    auto Xh = Pch<1>(mesh);


    auto thefms = fms( Xh, elements(mesh) );

    auto phio = Xh->element();
    phio = vf::project(Xh, elements(mesh), h() );
    phio +=vf::project(Xh, boundaryfaces(mesh), -idv(phio) - h()/100. );
    auto phi = thefms->operator()(phio);

    auto exp = exporter(_mesh=mesh, _name="disttowalls");
    exp->step(0)->add("phi", phi);
    exp->save();

}

int main( int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv, _about=makeAboutFSI());
    run();
    return 0;
}
