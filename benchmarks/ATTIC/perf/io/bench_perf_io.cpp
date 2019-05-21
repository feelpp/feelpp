#include <feel/feel.hpp>

int main(int argc, char**argv )
{

    double t0, t1, t2, t3, t4, t5;
    // MPI_Init ( &argc, &argv );
    t0 = MPI_Wtime();

    using namespace Feel;
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
        ( "mu", po::value<double>()->default_value( 1.0 ), "coeff" )
        ;

    Environment env( _argc=argc, _argv=argv,
                     _desc=laplacianoptions,
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM>>);

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();


    auto Vh = Pch<2>( mesh );
    auto u = Vh->element("u");
    auto mu = doption(_name="mu");
    auto f = expr( soption(_name="functions.f"), "f" );
    auto g = expr( soption(_name="functions.g"), "g" );
    auto v = Vh->element( g, "g" );

    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();

    //deleted calculations here
    
    //

    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "u", u );
    e->add( "g", v );

    MPI_Barrier(MPI_COMM_WORLD);
    t4 = MPI_Wtime();

    e->save();

    MPI_Barrier(MPI_COMM_WORLD);
    t5 = MPI_Wtime();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0)
    {
        const long dof = Vh->nDof();
        std::cout << "Laplacian program finished in " << (t5-t0)  << " seconds\n"
                  << "   degrees of freedom = " << dof << "\n" 
                  << "   init (?) time : " << (t1-t0) <<" seconds\n"
                  << "   loadMesh time : " << (t2-t1) <<" seconds\n"
                  << "   model loading time : " << (t3-t2) <<" seconds\n"
                  << "   exporter preparation time : " << (t4-t3) <<" seconds\n"
                  << "   writing time : " << (t5-t4) <<" seconds\n";
    }


    return 0;
    //# endmarker4 #
}
//! [global]

