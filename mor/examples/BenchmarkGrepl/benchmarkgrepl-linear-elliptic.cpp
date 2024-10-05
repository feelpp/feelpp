#include "benchmarkgrepl-linear-elliptic.hpp"

#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{
template<int Order>
BenchmarkGreplLinearElliptic<Order>::BenchmarkGreplLinearElliptic()
    :
    super_type( "BenchMarkGreplLinearElliptic" + std::to_string(Order) )
{
    this->setPluginName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_NAME) + fmt::format("P{}",Order) );
    this->setPluginLibName( BOOST_PP_STRINGIZE(FEELPP_MOR_PLUGIN_LIBNAME) );
}

template<int Order>
void BenchmarkGreplLinearElliptic<Order>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{

    M_mu = this->Dmu->element();

    std::shared_ptr<space_type> Xh;
    if ( !pT )
        pT.reset( new element_type );

    auto const& ptreeEim = ptree.get_child( "eim" );
    auto const& ptreeEimg = ptreeEim.get_child( "eim_g" );
    std::string dbnameEimg = ptreeEimg.template get<std::string>( "database-filename" );

    auto eim_g = eim( _model=eim_no_solve(super_type::shared_from_this()),
                      _element=*pT,
                      _space=Xh,
                      _parameter=M_mu,
                      _expr=1./sqrt( (Px()-cst_ref(M_mu(0)))*(Px()-cst_ref(M_mu(0))) + (Py()-cst_ref(M_mu(1)))*(Py()-cst_ref(M_mu(1))) ),
                      //_sampling=Pset,
                      _name="eim_g",
                      _filename=dbnameEimg,
                      _directory=dbDir );
    this->addEim( eim_g );
}

template<int Order>
void BenchmarkGreplLinearElliptic<Order>::initModel()
{

    std::string mshfile_name = soption(_name="mshfile");

    /*
     * First we create the mesh or load it if already exist
     */

    if( mshfile_name=="" )
    {
        double hsize=doption(_name="hsize");
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                               _desc = createGeo( hsize ) );

        //mesh = unitSquare( hsize );
    }
    else
    {
        mesh = loadMesh( _mesh=new mesh_type,
                         _filename=mshfile_name,
                         _update=MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
        // mesh = loadGMSHMesh( _mesh=new mesh_type,
        //                      _filename=option("mshfile").as<std::string>(),
        //                      _update=MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    }

    /*
     * The function space and some associate elements are then defined
     */
    auto Xh = space_type::New( mesh );
    this->setFunctionSpaces( Xh );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << Xh->nDof() << std::endl;
        std::cout << "Number of local dof " << Xh->nLocalDof() << std::endl;
    }

    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    this->Dmu->setDimension( 2 );
    auto mu_min = this->Dmu->element();
    mu_min <<  -1, -1;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << -0.01, -0.01;
    this->Dmu->setMax( mu_max );

    M_mu = this->Dmu->element();

    auto Pset = this->Dmu->sampling();
    //specify how many elements we take in each direction
    std::vector<size_type> N(2);
    int Ne = ioption(_name="trainset-eim-size");
    //40 elements in each direction
    N[0]=Ne; N[1]=Ne;

    //interpolation points are located on different proc
    //so we can't distribute parameters on different proc as in crb case
    //else for a given mu we are not able to evaluate g at a node wich
    //is not on the same proc than mu (so it leads to wrong results !)
    bool all_proc_same_sampling=true;
    std::string sampling_mode = "log-equidistribute";
    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();
    // std::string file_name = ( boost::format("eim_trainset_Ne%1%-proc%2%on%3%") % Ne %proc_number %total_proc).str();
    std::ifstream file ( supersamplingname );

    if( ! file )
    {
        //std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();
        //Pset->equidistributeProduct( N , all_proc_same_sampling , supersamplingname );
        //Pset->writeOnFile( supersamplingname );
        Pset->sample( N, sampling_mode, all_proc_same_sampling, supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto eim_g = eim( _model=eim_no_solve(super_type::shared_from_this()),
                      _element=*pT,
                      _space=Xh,
                      _parameter=M_mu,
                      _expr=1./sqrt( (Px()-cst_ref(M_mu(0)))*(Px()-cst_ref(M_mu(0))) + (Py()-cst_ref(M_mu(1)))*(Py()-cst_ref(M_mu(1))) ),
                      _sampling=Pset,
                      _name="eim_g" );
    this->addEim( eim_g );

#if 0
    if( Environment::worldComm().isMasterRank() )
    {
        bool error;
        std::cout<<" eim g mMax : "<<eim_g->mMax(error)<<" error : "<<error<<std::endl;
    }
#endif
    //checkEimExpansion();

    this->M_betaAqm.resize( 2 );
    this->M_betaAqm[0].resize( 1 );
    int eim_g_size = eim_g->mMax();
    this->M_betaAqm[1].resize( eim_g_size );
    this->M_betaFqm.resize( 2 );
    this->M_betaFqm[0].resize( 1 );
    this->M_betaFqm[0][0].resize( eim_g_size );
    this->M_betaFqm[1].resize( 1 );
    this->M_betaFqm[1][0].resize( 1 );


    assemble();

} // BenchmarkGreplLinearElliptic::init


template<int Order>
void BenchmarkGreplLinearElliptic<Order>::checkEimExpansion()
{
    auto Xh=this->Xh;
    auto Pset = this->Dmu->sampling();
    std::vector<size_type> N(2);
    int Ne = ioption(_name="trainset-eim-size");
    N[0]=Ne; N[1]=Ne;
    bool all_proc_same_sampling=true;
    std::string supersamplingname =(boost::format("PsetCheckEimExpansion-Ne%1%-generated-by-master-proc") %Ne ).str();
    std::ifstream file ( supersamplingname );
    if( ! file )
    {
        Pset->equidistributeProduct( N , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto eim_g = this->scalarContinuousEim()[0];

    //check that eim expansion of g is positive on each vertice
    int max = eim_g->mMax();
    auto e = exporter( _mesh=mesh );
    for( auto mu : *Pset )
    {

        if( Environment::worldComm().isMasterRank() )
            std::cout<<"check gM for mu = ["<< mu(0)<<" , "<<mu(1)<<"]"<<std::endl;

        auto exprg = 1./sqrt( (Px()-mu(0))*(Px()-mu(0)) + (Py()-mu(1))*(Py()-mu(1)) );
        auto g = vf::project( _space=Xh, _expr=exprg );

        for(int m=1; m<max; m++)
        {
            vectorN_type beta_g = eim_g->beta( mu , m );

            auto gM = expansion( eim_g->q(), beta_g , m);
            auto px=vf::project(_space=Xh, _expr=Px() );
            auto py=vf::project(_space=Xh, _expr=Py() );

            int size=Xh->nLocalDof();
            for(int v=0; v<size; v++)
            {
                double x = px(v);
                double y = py(v);
                if( gM(v) < -1e-13 )
                    std::cout<<"gM("<<x<<","<<y<<") = "<<gM(v)<<" ! - proc  "<<Environment::worldComm().globalRank()<<" but g("<<x<<","<<y<<") =  "<<g(v)<<" === m : "<<m<<std::endl;
                if( g(v) < -1e-13  )
                    std::cout<<"g("<<x<<","<<y<<") = "<<g(v)<<" donc la projection de g est negative"<<std::endl;
            }
#if 0

            for(int i=0; i<beta_g.size(); i++)
            {
                if( beta_g(i) < - 1e-14 && Environment::worldComm().isMasterRank() )
                {
                    std::cout<<"beta("<<i<<") is negative : "<<beta_g(i)<<"  === m : "<<m<<std::endl;
                }
            }
            //std::string name =( boost::format("GM%1%") %m ).str();
            //e->add( name , gM );
            //if( m == max-1 )
            //e->add( "exact" , g );

                //auto basis = eim_g->q(i);
                //double basisnorm = basis.l2Norm();
                //if( basisnorm < 1e-14 && Environment::worldComm().isMasterRank() )
                    //{
                    //std::cout<<"basis "<<i<<" norm negative : "<<basisnorm<<std::endl;
                    //exit(0);
                    //}
                //for(int v=0; v<size; v++)
                //{
                //    if( basis(i) < 1e-14 )
                //        std::cout<<"basis("<<v<<") = "<<basis(v)<<" ! - proc  "<<Environment::worldComm().globalRank()<<std::endl;
                //}
            }
#endif
        }//loop over m
        //e->save();
    }
}

template<int Order>
typename BenchmarkGreplLinearElliptic<Order>::monolithic_type
BenchmarkGreplLinearElliptic<Order>::computeMonolithicFormulation( parameter_type const& mu )
{
    auto Xh = this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma_dir = doption(_name="gamma");

    auto exprg = 1./sqrt( (Px()-mu(0))*(Px()-mu(0)) + (Py()-mu(1))*(Py()-mu(1)) );
    auto g = vf::project( _space=Xh, _expr=exprg );
    M_monoA = backend()->newMatrix( _test=Xh, _trial=Xh );
    M_monoF.resize(2);
    M_monoF[0] = backend()->newVector( Xh );
    M_monoF[1] = backend()->newVector( Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) = integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) + idt( u )*id( v )*idv(g) ) +
        integrate( _range=markedfaces( mesh, "boundaries" ), _expr=gamma_dir*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );
    form1( _test=Xh, _vector=M_monoF[0] ) = integrate( _range=elements(mesh) , _expr=id( v ) * idv(g) );
    form1( _test=Xh, _vector=M_monoF[1] ) = integrate( _range=elements(mesh) , _expr=id( v ) );

    return boost::make_tuple( M_monoA, M_monoF );

}


template<int Order>
void BenchmarkGreplLinearElliptic<Order>::assemble()
{
    auto Xh=this->Xh;
    auto u = Xh->element();
    auto v = Xh->element();

    double gamma_dir = doption(_name="gamma");
    auto eim_g = this->scalarContinuousEim()[0];

    this->M_Aqm.resize( 2 );
    this->M_Aqm[0].resize( 1 );
    this->M_Aqm[0][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range=markedfaces( mesh, "boundaries" ), _expr=gamma_dir*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );

    bool error;
    int M_g = eim_g->mMax(error);
    if( error ) M_g++;
    this->M_Aqm[1].resize( M_g );
    for(int m=0; m<M_g; m++)
    {
        this->M_Aqm[1][m] = backend()->newMatrix( _test=Xh, _trial=Xh );
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[1][m] ) =
            integrate( _range=elements( mesh ),
                       _expr=idt( u )* id( v ) * idv( eim_g->q(m) ) );
    }

    this->M_Fqm.resize( 2 );
    this->M_Fqm[0].resize( 1 );
    this->M_Fqm[1].resize( 1 );

    this->M_Fqm[0][0].resize(M_g);
    for(int m=0; m<M_g; m++)
    {
        this->M_Fqm[0][0][m] = backend()->newVector( Xh );
        form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) = integrate( _range=elements( mesh ),
                                                                     _expr=id( v ) * idv( eim_g->q(m) ) );
    }
    this->M_Fqm[1][0].resize(1);
    this->M_Fqm[1][0][0] = backend()->newVector( Xh );
    form1( _test=Xh, _vector=this->M_Fqm[1][0][0] ) = integrate( _range=elements( mesh ),
                                                                 _expr=id( v ) );


    //use computeLinearDecompositionA to provide innter product matrix
    //because this model use EIM
#if 0
    auto mu = refParameter();

    M = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) )+
        integrate( _range=markedfaces( mesh, "boundaries" ), _expr=gamma_dir*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u ) );
    this->addEnergyMatrix( M );

    vectorN_type beta_g = eim_g->beta( mu );
    for(int m=0; m<M_g; m++)
    {
        auto q = eim_g->q(m);
        q.scale( beta_g(m) );
        form2( _test=Xh, _trial=Xh, _matrix=M ) +=  integrate( _range=elements( mesh ), _expr= idt( u )*id( v ) * idv( q ) );
    }
#endif

}


template<int Order>
typename BenchmarkGreplLinearElliptic<Order>::vector_sparse_matrix
BenchmarkGreplLinearElliptic<Order>::computeLinearDecompositionA()
{
    auto Xh=this->Xh;
    auto muref = refParameter();
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma_dir = doption(_name="gamma");

    auto exprg = 1./sqrt( (Px()-muref(0))*(Px()-muref(0)) + (Py()-muref(1))* (Py()-muref(1)) );
    auto g = vf::project( _space=Xh, _expr=exprg );
    this->M_linearAqm.resize(1);
    this->M_linearAqm[0].resize(1);
    this->M_linearAqm[0][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] ) = integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) + idt( u )*id( v )*idv(g) ) +
        integrate( _range=markedfaces( mesh, "boundaries" ), _expr=gamma_dir*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );
    return this->M_linearAqm;

}

template<int Order>
typename BenchmarkGreplLinearElliptic<Order>::affine_decomposition_type
BenchmarkGreplLinearElliptic<Order>::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}


template<int Order>
double BenchmarkGreplLinearElliptic<Order>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve )
{

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    //solve the problem without affine decomposition
    //solution=this->solve(mu);

    double s=0;
    if ( output_index==0 )
    {
        for ( int q=0; q<1; q++ )
        {
            for ( int m=0; m<mMaxF(output_index,q); m++ )
            {
                s += this->M_betaFqm[output_index][q][m]*dot( *this->M_Fqm[output_index][q][m] , solution );
            }
        }
    }
    if( output_index==1 )
    {
        s = integrate( _range=elements( mesh ), _expr=idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}



template class BenchmarkGreplLinearElliptic<1>;
template class BenchmarkGreplLinearElliptic<2>;
template class BenchmarkGreplLinearElliptic<3>;

FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplLinearEllipticP1, BenchmarkGreplLinearElliptic<1>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,P1) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplLinearEllipticP2, BenchmarkGreplLinearElliptic<2>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,P2) )
FEELPP_CRB_PLUGIN_TEMPLATE( BenchmarkGreplLinearEllipticP3, BenchmarkGreplLinearElliptic<3>, BOOST_PP_CAT(FEELPP_MOR_PLUGIN_NAME,P3) )
}
