#include <boost/dll/alias.hpp> // for BOOST_DLL_ALIAS   
#include <feel/feelcrb/crbplugin.hpp>
#include <linearelasticity3d.hpp>

namespace Feel
{

po::options_description
makeLinearElasticity3dOptions()
{
    po::options_description solidoptions( "LinearElasticity3d options" );
    solidoptions.add_options()
        ( "gamma", po::value<double>()->default_value( 100 ), "penalisation term" )
        ( "poission-coeff", po::value<double>()->default_value( 0.3 ), "penalisation term" )
        ;
    return solidoptions;
}
AboutData
makeLinearElasticity3dAbout( std::string const& str )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "Linear Elasticity 3D Application",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2016 Feel++ Consortium" );
    return about;
}

LinearElasticity3d::LinearElasticity3d()
    :
    super_type( "LinearElasticity" )
{
    
}

void
LinearElasticity3d::initBetaQ()
{
    this->M_betaAq.resize( 4 );
    this->M_betaFq.resize( 2 );
    this->M_betaFq[0].resize( 1 );
    this->M_betaFq[1].resize( 1 );
}

LinearElasticity3d::super_type::betaq_type
LinearElasticity3d::computeBetaQ( parameter_type const& mu )
{
    //std::cout << "computeBetaQ start \n";
    if ( this->M_betaAq.empty() )
        this->initBetaQ();

    for (int k = 0;k<4;++k)
        this->M_betaAq[k] = mu(k);
    for (int k = 0;k<1;++k)
        this->M_betaFq[0][k] = mu( 4+k );
    this->M_betaFq[1][0] = 1.;
    //std::cout << "computeBetaQ finish \n";
    return boost::make_tuple( this->M_betaAq, this->M_betaFq );
}

void
LinearElasticity3d::initModel()
{

    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";

    typedef LinearElasticity3dConfig::mesh_type mesh_type;
    auto mesh = loadMesh( _mesh=new mesh_type,
                          _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES),
                          _savehdf5=0 );
    if( Environment::worldComm().isMasterRank() )
        std::cout << "mesh load with " << mesh->numGlobalElements() << " elements\n";

    this->setFunctionSpaces( LinearElasticity3dConfig::space_type::New( mesh ) );
    this->setSymbolicExpressionBuildDir("$repository/crb/linearelasticity3d/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";


    auto u = Xh->element();
    auto v = Xh->element();
    //auto mesh = Xh->mesh();

    double penaldir= doption(_name="gamma");
    //std::list<std::string> markersDispImposed = {"base1","base2","base3"};
    std::vector<std::pair<std::string,std::string> > markersBases = { { "mat-base1", "base1" },{ "mat-base2", "base2" },{ "mat-base3", "base3" } };

    /// [parameters]
    Dmu->setDimension( 5 );
    auto mu_min = Dmu->element();
    //mu_min << 1.e5,1.e5,1.e5,1.e5,1. ;
    mu_min << 1.e6,1.e6,1.e6,1.e6,1.e3 ;
    //mu_min << 1.e6,1.e2 ;
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    //mu_max << 1.e9,1.e9,1.e9,1.e9,5.e4;
    mu_max << 1.e8,1.e8,1.e8,1.e8,1.e4;
    //mu_max << 1.e6,5.e4;
    Dmu->setMax( mu_max );
    //double E_Ref = 1.e7;
    std::vector<double> E_Ref = { 1.1601e+06, 1.03232e+06, 1.42165e+06, 1.01663e+06 };

    //double E=1.e6;
    double nu=0.3;
    double lambda = nu/( (1+nu)*(1-2*nu) ); //E*nu/( (1+nu)*(1-2*nu) );
    double mu = 1./(2*(1+nu));// E/(2*(1+nu));
    auto Id = eye<mesh_type::nDim,mesh_type::nDim>();
    auto epsilon=sym(grad(u));
    auto epsilont=sym(gradt(u));
    auto sigma= lambda*trace(epsilon)*Id+2*mu*epsilon;
    auto sigmat= lambda*trace(epsilont)*Id+2*mu*epsilont;

    auto energy = backend()->newMatrix(_test=Xh,_trial=Xh);


    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate(_range=markedelements(mesh,"mat-central"),
                   _expr=inner(sigmat,grad(u)) );
    a0.matrixPtr()->close();
    this->addLhs( { a0 , "mu0" } );
    //energy->addMatrix(1.,a0.matrixPtr() );
    energy->addMatrix(E_Ref[0],a0.matrixPtr() );

    for (int k=0;k<3;++k )
    {
        std::string const& markerVol = markersBases[k].first;
        std::string const& markerSurf = markersBases[k].second;
        auto a1 = form2( _trial=Xh, _test=Xh);
        a1 = integrate(_range=markedelements(mesh,markerVol),
                       _expr=inner(sigmat,grad(u)) );
        a1 += integrate(_range=markedfaces(mesh,markerSurf),
                        _expr=-inner(sigmat*N(),id(u)) + inner(-sigma*N()+std::max(2*mu,lambda)*penaldir*id(u)/hFace(),idt(u)) );
        a1.matrixPtr()->close();
        this->addLhs( { a1 , (boost::format("mu%1%")%(k+1)).str() } );
        //energy->addMatrix(1.,a1.matrixPtr() );
        energy->addMatrix( E_Ref[k+1],a1.matrixPtr() );
    }

    //int cptMark=4;
    auto f1 = form1( _test=Xh );
    f1 = integrate( _range=markedfaces( mesh,"support-top" ), _expr=/*-1.e4**/-inner(id(u),N()) );
    f1.vectorPtr()->close();
    //this->addRhs( { f1, (boost::format("mu%1%")%cptMark++).str() } );
    //this->addRhs( { f1, "1" } );
    this->addRhs( { f1, "mu4" } );



    /// [energy]
    //a0.matrixPtr()->symmetricPart( energy );
    //energy->addMatrix(1.,a0.matrixPtr() );
#if 1
    energy->close();
    this->addEnergyMatrix( energy );
#else
    auto energy2 = form2( _trial=Xh, _test=Xh);
    energy2 = integrate(_range=elements(mesh),
                        _expr=E_Ref[0]*/*mu**/inner(gradt(u),grad(u)) );
    this->addEnergyMatrix( energy2 );

#endif

    /// [output]
    auto out1 = form1( _test=Xh );
#if 0
    double meas = integrate(_range=elements(mesh),_expr=cst(1.)).evaluate()(0,0);
    out1 = integrate( _range=elements( mesh ), _expr=div( u )/cst(meas)) ;
#else
    double meas = integrate(_range=elements(mesh),_expr=cst(1.)).evaluate()(0,0);
    out1 = integrate( _range=elements( mesh ), _expr=(id( u )(1,0))/cst(meas)) ;
    out1.vectorPtr()->close();
#endif
    this->addOutput( { out1, "1" } );

#if 1
    auto mode1 = Xh->element( oneX() );
    auto mode2 = Xh->element( oneY() );
    auto mode3 = Xh->element( oneZ() );
    auto mode4 = Xh->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = Xh->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = Xh->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
    std::shared_ptr<NullSpace<value_type> > myNearNullSpace( new NullSpace<value_type>(backend( _name="backend-primal"),userNullSpace) );
#else
    auto userNullSpace = nullspace_ptr( Xh, oneX(),oneY(),oneZ(),vec(Py(),-Px(),cst(0.)),vec(-Pz(),cst(0.),Px()),vec(cst(0.),Pz(),-Py()) );
#endif
    backend( _name="backend-primal")->attachNearNullSpace( myNearNullSpace );
    backend( _name="backend-dual")->attachNearNullSpace( myNearNullSpace );
    backend( _name="backend-l2")->attachNearNullSpace( myNearNullSpace );

}

double
LinearElasticity3d::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    //CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;

    if ( output_index == 1 )
    {
        double meas = integrate(_range=elements(mesh),_expr=cst(1.)).evaluate()(0,0);
#if 0
        output = integrate(_range=elements(mesh),_expr=divv(u)/cst(meas)).evaluate()(0,0);
#else
        output = integrate(_range=elements(mesh),_expr=(idv(u)(1,0))/cst(meas)).evaluate()(0,0);
#endif
        std::cout << " LinearElasticity3d::output " << output << "\n";
    }
    else
        throw std::logic_error( "[LinearElasticity3d::output] error with output_index : only 1 " );

    return output;
}

FEELPP_CRB_PLUGIN( LinearElasticity3d, linearelasticity3d )
} // namespace Feel
