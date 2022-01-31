#include "benchmarkgrepl-nonlinear-elliptic.hpp"

#include <feel/feelmor/crbplugin.hpp>

namespace Feel
{

template<int Order, int Dim>
BenchmarkGreplNonlinearElliptic<Order,Dim>::BenchmarkGreplNonlinearElliptic(std::string const& prefix)
    :
    super_type( "BenchMarkGreplNonlinearElliptic" + std::to_string(Order) + "_" + std::to_string(Dim) + "D" ),
    M_use_newton( boption(_name="crb.use-newton") ),
    M_useSerErrorEstimation( boption(_name="ser.error-estimation") ),
    M_use_deim( boption("use-deim") )
{}


template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    M_mu = this->Dmu->element();

    std::shared_ptr<space_type_eimg> Xh_eimg;
    if ( !pT )
        pT.reset( new element_type );

    auto const& ptreeEim = ptree.get_child( "eim" );
    auto const& ptreeEimg = ptreeEim.get_child( "eim_g" );
    std::string dbnameEimg = ptreeEimg.template get<std::string>( "database-filename" );

    auto eim_g = eim( _model=std::dynamic_pointer_cast< BenchmarkGreplNonlinearElliptic<Order,Dim> >( this->shared_from_this() ),
                      _element=*pT,
                      _space=Xh_eimg,
                      _parameter=M_mu,
                      _expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*_e1 ) - 1 ),
                      //_sampling=Pset,
                      _name="eim_g",
                      _filename=dbnameEimg,
                      _directory=dbDir );
    this->addEim( eim_g );
}
template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::initModel()
{

    if ( !this->Xh )
    {
        /*
         * First we create the mesh or load it if already exist
         */
        static mesh_ptrtype static_mesh;
        if ( !static_mesh )
        {
            std::string mshfile_name = soption(_name="mshfile");
            if( mshfile_name=="" )
            {
                double hsize=doption(_name="hsize");
                mesh = createGMSHMesh( _mesh=new mesh_type,
                                       _desc=domain( _name = "benchmarkgrepl",
                                                     _shape = "hypercube",
                                                     _dim = Dim,
                                                     _h=hsize,
                                                     _xmin=0,_xmax=1,
                                                     _ymin=0,_ymax=1 ) );
            }
            else
            {
                mesh = loadMesh( _mesh=new mesh_type,
                                 _filename=soption("mshfile") );
            }
            static_mesh = mesh;
        }
        else
            mesh = static_mesh;
        /*
         * The function space and some associate elements are then defined
         */
        this->setFunctionSpaces( space_type::New( mesh ) );
        // allocate an element of Xh
        pT = element_ptrtype( new element_type( this->Xh ) );
    }

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << std::endl;
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << std::endl;
    }


    this->Dmu->setDimension( 2 );
    auto mu_min = this->Dmu->element();
    mu_min <<  0.01, 0.01;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    mu_max << 10, 10;
    this->Dmu->setMax( mu_max );

    M_mu = this->Dmu->element();

    auto Pset = this->Dmu->sampling();
    int Ne = ioption(_name="trainset-eim-size");
    std::string supersamplingname =(boost::format("DmuEim-Ne%1%-generated-by-master-proc") %Ne ).str();

    std::ifstream file ( supersamplingname );
    if( !file )
    {
        //specify how many elements we take in each direction
        std::vector<size_type> N(2);
        //40 elements in each direction
        N[0]=Ne; N[1]=Ne;

        //interpolation points are located on different proc
        //so we can't distribute parameters on different proc as in crb case
        //else for a given mu we are not able to evaluate g at a node wich
        //is not on the same proc than mu (so it leads to wrong results !)
        bool all_proc_same_sampling=true;
        std::string samplingMode = "equidistribute";
        Pset->sample( N, samplingMode, all_proc_same_sampling , supersamplingname );
        //Pset->equidistributeProduct( N , all_proc_same_sampling , supersamplingname );
        //Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    int M = 0;
    if ( M_use_deim )
    {
        auto d = deim( _model=std::dynamic_pointer_cast<self_type>(this->shared_from_this()),
                       _sampling=Pset );
        this->addDeim( d );
        this->deim()->run();
        M = this->deim()->size();
    }
    else
    {
        space_ptrtype_eimg Xh_eimg = space_type_eimg::New( mesh );
        auto eim_g = eim( _model=std::dynamic_pointer_cast<self_type >( this->shared_from_this() ),
                          _element=*pT,
                          _space=Xh_eimg,
                          _parameter=M_mu,
                          //_expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*idv(*pT) ) - 1 ),
                          _expr=( cst_ref(M_mu(0))/cst_ref(M_mu(1)) )*( exp( cst_ref(M_mu(1))*_e1 ) - 1 ),
                          _sampling=Pset,
                          _name="eim_g" );

        this->addEim( eim_g );
        M = this->scalarContinuousEim()[0]->mMax();
    }

    this->initBetaQm( M );

    //resize data structure
    if( M_use_newton )
    {
        this->M_Jqm.resize( 3 );
        this->M_Jqm[0].resize( 1 );
        this->M_Jqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        this->M_Jqm[1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        }
        this->M_Jqm[2].resize( 1 );
        this->M_Jqm[2][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Rqm.resize( 2 );
        this->M_Rqm[0].resize( 3 );
        this->M_Rqm[1].resize( 1 );
        this->M_Rqm[0][0].resize(1);
        this->M_Rqm[0][0][0] = backend()->newVector( this->Xh );
        this->M_Rqm[0][1].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        }
        this->M_Rqm[0][2].resize(1);
        this->M_Rqm[0][2][0] = backend()->newVector( this->Xh );
        this->M_Rqm[1][0].resize( 1 );
        this->M_Rqm[1][0][0]= backend()->newVector( this->Xh );
    }

    if( !M_use_newton || M_useSerErrorEstimation )
    {
        this->M_Aqm.resize( 1 );
        this->M_Aqm[0].resize( 1 );
        this->M_Aqm[0][0] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );

        this->M_Fqm.resize( 2 );
        this->M_Fqm[0].resize( 2 );
        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
        }
        this->M_Fqm[0][1].resize(1);
        this->M_Fqm[0][1][0]=backend()->newVector( this->Xh );
        this->M_Fqm[1].resize(1);
        this->M_Fqm[1][0].resize(1);
        this->M_Fqm[1][0][0]=backend()->newVector( this->Xh );
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = this->Xh->elementPtr();

    assemble();

} // BenchmarkGreplNonlinearElliptic::init

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::assembleJacobianWithAffineDecomposition( std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto u = Xh->element(); //trial
    double gamma = doption(_name="gamma");

    int M;
    if ( M_use_deim )
        M = this->deim()->size();
    else
        M = this->scalarContinuousEim()[0]->mMax();

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=Jqm[0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idt(u)*id(v)/hFace()
                   - (gradt(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idt(u) );

    Jqm[0][0]->close();

    Jqm[1].resize(M);
    for(int m=0; m<M; m++)
    {
        this->M_Jqm[1][m] = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        if ( M_use_deim )
        {
            v = *(this->deim()->q(m));
            form2( _test=Xh, _trial=Xh, _matrix=Jqm[1][m] ) =
                integrate( _range = elements(mesh), _expr = idv(v)*idt(u)*id(v) );
        }
        else
        {
            auto eim_g = this->scalarContinuousEim()[0];
            form2( _test=Xh, _trial=Xh, _matrix=Jqm[1][m] ) =
                integrate( _range = elements(mesh), _expr = idv(eim_g->q(m))*idt(u)*id(v) );

        }
        Jqm[1][m]->close();
    }

    form2( _test=Xh, _trial=Xh, _matrix=Jqm[2][0] ) =
        integrate( _range= elements( mesh ), _expr = idt(u)*id(v) );
    Jqm[2][0]->close();

}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto Xh = this->Xh;
    auto v = Xh->element(); //test
    auto u = Xh->element();
    double gamma = doption(_name="gamma");

    int M;
    if ( M_use_deim )
        M = this->deim()->size();
    else
        M = this->scalarContinuousEim()[0]->mMax();

    u = *this->M_InitialGuess[0][0];

    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );
    Rqm[0][0][0]->close();


    Rqm[0][1].resize(M);
    for(int m=0; m<M; m++)
    {
        this->M_Rqm[0][1][m] = backend()->newVector( this->Xh );
        if ( M_use_deim )
        {
            v = *(this->deim()->q(m));
            form1( _test=Xh, _vector=Rqm[0][1][m] ) =
                integrate( _range= elements( mesh ), _expr=( idv(v)*id(v) ) );
        }
        else
        {
            auto eim_g = this->scalarContinuousEim()[0];
            form1( _test=Xh, _vector=Rqm[0][1][m] ) =
                integrate( _range= elements( mesh ), _expr=( idv(eim_g->q(m))*id(v) ) );
        }
        Rqm[0][1][m]->close();
    }

    if( Dim == 2 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=Rqm[0][2][0] ) =
            integrate( _range= elements( mesh ), _expr=sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=Rqm[1][0][0] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    Rqm[0][2][0]->close();
    Rqm[1][0][0]->close();
}

template <int Order, int Dim>
bool
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm)
{
    auto Xh = this->Xh;
    auto u = Xh->element();
    auto v = Xh->element(); //test
    double gamma = doption(_name="gamma");

    u = X;
    Rqm[0][0][0] = backend()->newVector( this->Xh );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );
    form1( _test=Xh, _vector=Rqm[0][0][0] ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    Rqm[0][0][0]->close();
    return true;
}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateJacobianMonolithic( vector_ptrtype const& X,
                                                                  sparse_matrix_ptrtype & J,
                                                                  parameter_type const& mu )
{
    //Here we don't use EIM approximation and so no affine decomposition
    auto Xh = this->Xh;
    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    double gamma = doption(_name="gamma");
    auto g = exp( mu(1)*idv(u) );

    J->zero();
    form2( _test=Xh, _trial=Xh, _matrix=J ) =
        integrate( _range= elements( mesh ),_expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = boundaryfaces( mesh ),
                                                          _expr = gamma*idt(u)*id(v)/hFace()
                                                          - (gradt(u)*vf::N())*id(v)
                                                          - (grad(v)*vf::N())*idt(u) );

    form2( _test=Xh, _trial=Xh, _matrix=J ) +=
        integrate( _range = elements(mesh), _expr = mu(0)*exp( mu(1)*idv(u) )*idt(u)*id(v) );

    J->close();

}

template <int Order, int Dim>
void
BenchmarkGreplNonlinearElliptic<Order,Dim>::updateResidualMonolithic(vector_ptrtype const& X,
                                                                 vector_ptrtype & R,
                                                                 parameter_type const& mu )
{
    auto Xh = this->Xh;
    double gamma = doption(_name="gamma");

    auto u = Xh->element();
    auto v = Xh->element(); //test
    u=*X;
    auto g = exp( mu(1)*idv(u) );

    R->zero();
    form1( _test=Xh, _vector=R ) =
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v)
                   - (grad(v)*vf::N())*idv(u) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= mu(0)/mu(1)*( exp( mu(1)*idv(u) ) * id(v) ) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr= -mu(0)/mu(1)*( id(v)) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if( Dim == 3 )
        form1( _test=Xh, _vector=R ) +=
            integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    R->close();

}

template<int Order, int Dim>
void BenchmarkGreplNonlinearElliptic<Order,Dim>::assemble()
{
    auto Xh = this->Xh;
    auto u = Xh->element();

    if( M_use_newton )
    {
        assembleJacobianWithAffineDecomposition( this->M_Jqm );
        assembleResidualWithAffineDecomposition( this->M_Rqm );
    }

    if( !M_use_newton || M_useSerErrorEstimation )
    {
        auto v = Xh->element();
        double gamma = doption(_name="gamma");
        int M;
        if ( M_use_deim )
            M = this->deim()->size();
        else
            M = this->scalarContinuousEim()[0]->mMax();

        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) =
            integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
        form2( _test=Xh, _trial=Xh, _matrix=this->M_Aqm[0][0] ) +=
            integrate( _range = boundaryfaces( mesh ),
                       _expr = gamma*idt(u)*id(v)/hFace()
                       - (gradt(u)*vf::N())*id(v)
                       - (grad(v)*vf::N())*idt(u) );

        this->M_Fqm[0][0].resize( M );
        for(int m=0; m<M; m++)
        {
            this->M_Fqm[0][0][m] = backend()->newVector( this->Xh );
            if ( M_use_deim )
            {
                v = *(this->deim()->q(m));
                form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) =
                    integrate( _range= elements( mesh ), _expr=( idv(v)*id(v) ) );
            }
            else
            {
                auto eim_g = this->scalarContinuousEim()[0];
                form1( _test=Xh, _vector=this->M_Fqm[0][0][m] ) =
                    integrate( _range= elements( mesh ), _expr=( idv(eim_g->q(m))*id(v) ) );
            }
            this->M_Fqm[0][0][m]->close();
        }

        if( Dim == 2 )
            form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
                integrate( _range= elements( mesh ),_expr=sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
        else if ( Dim == 3 )
            form1( _test=Xh, _vector=this->M_Fqm[0][1][0] ) =
                integrate( _range= elements( mesh ),_expr=sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );
        this->M_Fqm[0][1][0]->close();

        form1( _test=Xh, _vector=this->M_Fqm[1][0][0] ) =
            integrate( _range= elements( mesh ), _expr=id(v) );
        this->M_Fqm[1][0][0]->close();
    }//no newton

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::monolithic_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution)
{
    auto Xh = this->Xh;
    auto u = Xh->element();
    auto v = Xh->element();
    double gamma = doption(_name="gamma");

    M_monoA = backend()->newMatrix( _test=Xh, _trial=Xh );
    M_monoF.resize(2);
    M_monoF[0] = backend()->newVector( this->functionSpace() );
    M_monoF[1] = backend()->newVector( Xh );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) =
        integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=M_monoA ) += integrate( _range = boundaryfaces( mesh ),
                                                                _expr = gamma*idt(u)*id(v)/hFace()
                                                                - (gradt(u)*vf::N())*id(v)
                                                                - (grad(v)*vf::N())*idt(u) );

    if( Dim == 2 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );
    else if ( Dim == 3 )
        form1( _test=Xh, _vector=M_monoF[0] ) =
            integrate( _range= elements( mesh ), _expr=-mu(0)/mu(1)*( exp( mu(1)*idv(solution) ) - 1 )*id(v) + 100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * id(v) );

    form1( _test=Xh, _vector=M_monoF[1] ) =
        integrate( _range= elements( mesh ), _expr=id(v) );

    M_monoF[0]->close();
    M_monoF[1]->close();
    M_monoA->close();

    return boost::make_tuple( M_monoA, M_monoF );

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::element_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::solve( parameter_type const& mu )
{
    auto Xh=this->Xh;
    auto backendMonolithic = backend(_rebuild=true);
    sparse_matrix_ptrtype J = backendMonolithic->newMatrix( _test=Xh, _trial=Xh);
    vector_ptrtype R = backendMonolithic->newVector( Xh );
    backendMonolithic->nlSolver()->jacobian = boost::bind( &self_type::updateJacobianMonolithic,
                                                   boost::ref( *this ), _1, _2, mu );
    backendMonolithic->nlSolver()->residual = boost::bind( &self_type::updateResidualMonolithic,
                                                   boost::ref( *this ), _1, _2, mu );

    auto solution = Xh->element();
    backendMonolithic->nlSolve(_jacobian=J, _solution=solution, _residual=R);

    return solution;

}


template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::vector_sparse_matrix
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeLinearDecompositionA()
{
    auto Xh=this->Xh;
    auto u=Xh->element();
    auto v=Xh->element();
    double gamma = doption(_name="gamma");

    this->M_linearAqm.resize(1);
    this->M_linearAqm[0].resize(1);
    this->M_linearAqm[0][0] = backend()->newMatrix( _test=Xh, _trial=Xh );
    form2( _test=Xh, _trial=Xh, _matrix=this->M_linearAqm[0][0] ) =
        integrate( _range=elements( mesh ), _expr=gradt( u )*trans( grad( v ) ) ) +
        integrate( _range=boundaryfaces( mesh ), _expr=gamma*idt( u )*id( v )/h()
                   -gradt( u )*vf::N()*id( v )
                   -grad( v )*vf::N()*idt( u )
                   );

    this->M_linearAqm[0][0]->close();
    return this->M_linearAqm;

}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::affine_decomposition_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeAffineDecomposition()
{
    if( M_use_newton )
        return boost::make_tuple( this->M_Jqm, this->M_Rqm );
    else
        return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order, int Dim>
typename BenchmarkGreplNonlinearElliptic<Order,Dim>::affine_decomposition_type
BenchmarkGreplNonlinearElliptic<Order,Dim>::computePicardAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm );
}

template<int Order, int Dim>
std::vector< std::vector< typename BenchmarkGreplNonlinearElliptic<Order,Dim>::element_ptrtype > >
BenchmarkGreplNonlinearElliptic<Order,Dim>::computeInitialGuessAffineDecomposition( )
{
    return M_InitialGuess;
}


template<int Order, int Dim>
double BenchmarkGreplNonlinearElliptic<Order,Dim>::output( int output_index, parameter_type const& mu, element_type &solution, bool need_to_solve )
{
    //CHECK( need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;
    if ( output_index==0 )
    {
        if( Dim == 2 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * idv( solution ) ).evaluate()(0,0);
        else if ( Dim == 3 )
            s = integrate( _range= elements( mesh ), _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py())*sin(2*M_PI*Pz()) * idv( solution ) ).evaluate()(0,0);
    }
    if( output_index==1 )
    {
        s = integrate( _range=elements( mesh ), _expr=idv( solution ) ).evaluate()( 0,0 );
    }

    return s ;
}


}
