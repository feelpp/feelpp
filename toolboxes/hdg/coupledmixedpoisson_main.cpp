#include <feel/feelmodels/hdg/coupling.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "mixed-poisson-model" ,
                     "mixed-poisson-model" ,
                     "0.1",
                     "Mixed-Poisson-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016-2019 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Romain Hild", "developer", "", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;
}



// template<int nDim, int OrderT>
void
runApplicationCoupledMixedPoisson()
{
    using namespace Feel;

    const int nDim = 3;
    const int nOrder = 2;

    typedef FeelModels::CoupledMixedPoisson<nDim,nOrder> cmp_type;

    auto CMP = cmp_type::New();
    auto mesh = loadMesh( _mesh=new cmp_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<nOrder>(mesh), _imageSpace=Pdh<nOrder>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<nOrder>(mesh), _imageSpace=Pdhv<nOrder>(mesh) ) ) Idhv;
    if ( soption( "gmsh.submesh" ).empty() )
        CMP -> init(mesh);
    else
    {
        Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
		auto cmesh = createSubmesh( _mesh=mesh, _range=markedelements(mesh,soption("gmsh.submesh")) );
        Idh = IPtr( _domainSpace=Pdh<nOrder>(cmesh), _imageSpace=Pdh<nOrder>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhv<nOrder>(cmesh), _imageSpace=Pdhv<nOrder>(mesh) );
        CMP -> init( cmesh, mesh );
    }

    CMP -> run ( Idh, Idhv);

}


int main(int argc, char *argv[])
{
    using namespace Feel;

    po::options_description mpoptions( "hdg.poisson options" );
    mpoptions.add( FeelModels::makeCoupledMixedPoissonOptions("","hdg.poisson") );
    mpoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
        ;

    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=mpoptions,
                           _desc_lib=FeelModels::makeCoupledMixedPoissonLibOptions("","hdg.poisson").add(feel_options())
                           );

    /*
    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);
#if FEELPP_INSTANTIATION_ORDER_MAX >= 3
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ),
                                             hana::make_tuple("P3", hana::int_c<3> ) );
#elif FEELPP_INSTANTIATION_ORDER_MAX >= 2
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ),
                                             hana::make_tuple("P2", hana::int_c<2> ) );
#else
    auto discretizationt = hana::make_tuple( hana::make_tuple("P1", hana::int_c<1> ) );
#endif

    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension]( auto const& d )
                                                                                         {
                                                                                             constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                                                                                             std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                                                                                             constexpr int _torder = std::decay_t<decltype(hana::at_c<1>( hana::at_c<1>(d) ))>::value;
                                                                                             if ( dimension == _dim && discretization == _discretization )
                                                                                                 runApplicationCoupledMixedPoisson<_dim,_torder>();
                                                                                         } );
    */
    
    runApplicationCoupledMixedPoisson();

    return 0;
}
