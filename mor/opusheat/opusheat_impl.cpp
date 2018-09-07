
#include <opusheat.hpp>


namespace Feel
{

po::options_description
makeOpusHeatOptions()
{
    po::options_description heatequationoptions( "OpusHeat options" );
    heatequationoptions.add_options()
    ( "qsource", Feel::po::value<double>()->default_value( 1.e6 ), "qsource" )
    ( "vinconv", Feel::po::value<double>()->default_value( 1 ), "qsource" )
    ( "mu-max", Feel::po::value<double>()->default_value( 1.1 ), "qsource" )
    ( "do-not-use-operators-free", Feel::po::value<bool>()->default_value( true ), "never use operators free if true" )
    ;
    return heatequationoptions.add( bdf_options( "opusheat" ) ).add( ts_options( "opusheat" ) );
}
AboutData
makeOpusHeatAbout( std::string const& str )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "opusheat",
                           Feel::AboutData::License_GPL,
                           "Copyright (C) 2009-2014 Feel++ Consortium");
    return about;
}

template<bool IsStationary>
OpusHeat<IsStationary>::OpusHeat()
    :
    super_type((IsStationary)?"OpusHeat_stationary":"OpusHeat")
{}

template<bool IsStationary>
void
OpusHeat<IsStationary>::initDataStructureAffineDecomposition()
{
    M_Qa=5;
    M_Qm=1;
    M_Nl=2;
    M_Ql.resize( 2 );
    M_Ql[0]=4;
    M_Ql[1]=1;

    this->M_betaAq.resize( M_Qa );
    if ( super_type::is_time_dependent )
        this->M_betaMq.resize( M_Qm );
    this->M_betaFq.resize( M_Nl );
    for(int i=0; i<M_Nl; i++)
    {
        int ql=M_Ql[i];
        this->M_betaFq[i].resize( ql );
    }

    this->M_Aq.resize( M_Qa );
    if ( super_type::is_time_dependent )
        this->M_Mq.resize( M_Qm );
    this->M_Fq.resize( M_Nl );
    for(int l=0; l<M_Nl; l++)
    {
        this->M_Fq[l].resize( M_Ql[l] );
    }

}

template<bool IsStationary>
void
OpusHeat<IsStationary>::initializationField( element_ptrtype& initial_field,parameter_type const& mu )
{
    if ( super_type::is_time_dependent )
        initial_field->setConstant( 300. );
}

template<bool IsStationary>
void
OpusHeat<IsStationary>::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    auto ptreeSpecificityOfModel = ptree.get_child_optional( "specifity-of-model" );
    CHECK( ptreeSpecificityOfModel ) << "invalid ptree : section specifity-of-model is missing";
    M_measureMarkedSurface["IC2"] = ptreeSpecificityOfModel->template get<double>( "surface-measure-IC2" );
    //std::cout << "surface loaded " << M_surface << "\n";
}
template<bool IsStationary>
void
OpusHeat<IsStationary>::updateSpecificityModel( boost::property_tree::ptree & ptree ) const
{
    ptree.add( "surface-measure-IC2", M_measureMarkedSurface.find("IC2")->second );
}


template<bool IsStationary>
void
OpusHeat<IsStationary>::initModel()
{

    auto mesh = loadMesh( _mesh=new typename OpusHeatConfig<IsStationary>::mesh_type,
                          //_update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES),
                          _savehdf5=0 );
    this->setFunctionSpaces( space_type::New( mesh ) );

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of local dof " << this->Xh->nLocalDof() << "\n";
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }

    for ( std::string const& marker : std::vector<std::string>({"AIR","PCB","IC1","IC2"}) )
        M_measureMarkedSurface[marker] = measure(_range=markedelements(mesh,marker) );

    if ( super_type::is_time_dependent )
        M_bdf = bdf( _space=this->Xh, _vm=Environment::vm(), _name="opusheat" , _prefix="opusheat" );

    this->initDataStructureAffineDecomposition();

    this->Dmu->setDimension( 6 );
    auto mu_min = this->Dmu->element();
    mu_min << 1,1,0.1,1e4,1e4,1;
    this->Dmu->setMin( mu_min );
    auto mu_max = this->Dmu->element();
    //mu_max << 1.1;
    mu_max << 3,3, 5,1e6,1e6,30;
    //mu_max << doption(_name="mu-max");
    this->Dmu->setMax( mu_max );

    assembleData();

} // OpusHeat::init

template<bool IsStationary>
void
OpusHeat<IsStationary>::updateBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent )
{
    if ( this->M_betaAq.empty() )
        this->initDataStructureAffineDecomposition();

    double kappaIC1 = mu( 0 );
    double kappaIC2 = mu( 1 );
    double hcoeff = mu( 2 );
    double qIC1 = mu( 3 );
    double qIC2 = mu( 4 );
    double vinconv = mu( 5 );
    //if( !only_terms_time_dependent )
    {
        this->M_betaAq[0] = 1 ;
        this->M_betaAq[1] = kappaIC1 ;
        this->M_betaAq[2] = kappaIC2 ;
        this->M_betaAq[3] = vinconv;//1;//uconvRate  ;
        this->M_betaAq[4] = hcoeff;
        if ( super_type::is_time_dependent )
            this->M_betaMq[0] = 1;
    }
    this->M_betaFq[0][0] = 1;
    this->M_betaFq[0][1] = hcoeff;
    this->M_betaFq[0][2] = qIC1;
    this->M_betaFq[0][3] = qIC2;

    this->M_betaFq[1][0] = 1./M_measureMarkedSurface["IC2"];
}



template<bool IsStationary>
typename OpusHeat<IsStationary>::super_type::betaq_type
OpusHeat<IsStationary>::computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent )
{
    return computeBetaQ_impl<super_type>( mu, time, only_terms_time_dependent );
}
template<bool IsStationary>
typename OpusHeat<IsStationary>::super_type::betaq_type
OpusHeat<IsStationary>::computeBetaQ( parameter_type const& mu )
{
    return computeBetaQ_impl<super_type>( mu );
}
template<bool IsStationary>
void
OpusHeat<IsStationary>::assembleData()
{
    auto mesh = this->Xh->mesh();
    auto u = this->Xh->element();
    auto v = this->Xh->element();

    std::map<std::string,double> conductivity;
    conductivity["AIR"]=3e-2;
    conductivity["PCB"]=0.2;

    std::map<std::string,double> rhoC;
    rhoC["AIR"]=1100;//1;//1100;
    rhoC["PCB"]=2e6;//1;//2e6;
    rhoC["IC1"]=1.4e6;//1;//1.4e6;
    rhoC["IC2"]=1.4e6;//1;//1.4e6;

    double vinconv = 1.;//doption(_name="vinconv");
    double Qsource = doption(_name="qsource");//1.e6;

    if( true ) //boption("do-not-use-operators-free") )
    {
        for (int k = 0 ; k<this->M_Aq.size() ; ++k)
            this->M_Aq[k] = backend()->newMatrix( this->Xh, this->Xh );
        if ( super_type::is_time_dependent )
            this->M_Mq[0] = backend()->newMatrix( this->Xh, this->Xh );
        for (int k = 0 ; k<this->M_Fq[0].size() ; ++k)
            this->M_Fq[0][k] = backend()->newVector( this->Xh );
        this->M_Fq[1][0] = backend()->newVector( this->Xh );

        for ( std::string const& marker : std::vector<std::string>({"AIR","PCB"}) )
            form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
                integrate( _range= markedelements( mesh,marker ),
                           _expr= conductivity[marker]*gradt( u )*trans( grad( v ) ) );

        // weak Dirichlet condition (lhs)
        double penaldir=10.;
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
            integrate( _range=markedfaces( mesh, std::list<std::string>({"Gamma_4_AIR1","Gamma_4_AIR4"} ) ),
                       _expr= conductivity["AIR"]*( -(gradt(u)*id(v)+grad(v)*idt(u))*N() + penaldir*idt(u)*id(v)/hFace() ) );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0]) +=
            integrate( _range=markedfaces( mesh, "Gamma_4_PCB" ),
                       _expr= conductivity["PCB"]*( -(gradt(u)*id(v)+grad(v)*idt(u))*N() + penaldir*idt(u)*id(v)/hFace() ) );

        // parameter conductivity
        for ( std::string const& marker : std::vector<std::string>({"IC1"}) )
            form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[1])
                += integrate( _range= markedelements( mesh,marker ),
                              _expr= gradt( u )*trans( grad( v ) ) );
        // parameter conductivity
        for ( std::string const& marker : std::vector<std::string>({"IC2"}) )
            form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[2])
                += integrate( _range= markedelements( mesh,marker ),
                              _expr= gradt( u )*trans( grad( v ) ) );

        auto uconv = this->Xh->element( vinconv*expr("-(x-0.008)*(x-0.054)*(x>0.008):x" ) );
        auto uconvexpr = idv(uconv);
        //auto uconvexpr = vinconv*expr("-(x-0.008)*(x-0.054)*(x>0.008):x" );
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[3])
            += integrate( _range= markedelements( mesh,"AIR" ),
                          _expr= rhoC["AIR"]*(gradt(u)*uconvexpr*oneY())*id( v ) );

        // Robin condition in lhs (parameter coeff h)
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[4])
            += integrate( _range= markedfaces( mesh,"Gamma_2" ),
                          _expr= idt(u)*id(v) );

#if 0
        auto vDirExpr = rhoC["AIR"]*uconvexpr*oneY();
        //auto coeffStab = cst(1.);//vf::h()/( 2*norm2(  vDirExpr ) );
        auto chiConv = chi( norm2(  vDirExpr ) > 1e-6 );
        auto coeffStab = chiConv*vf::h()/( 2*(chiConv*norm2( vDirExpr)+(1-chiConv) ) );
        auto residualStabForm2 = gradt(u)*vDirExpr;
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[0/*2*/])
            += integrate(_range=markedelements(mesh,"AIR"),
                         _expr= coeffStab*(/*theta**/grad(v)*vDirExpr)*residualStabForm2 );
#endif
#if 0
        auto vDirExpr = rhoC["AIR"]*uconvexpr*oneY();
        //auto coeffStab = cst(1.);//vf::h()/( 2*norm2(  vDirExpr ) );
        auto chiConv = chi( norm2(  vDirExpr ) > 1e-6 );
        auto coeffStab = chiConv*vf::h()/( 2*(chiConv*norm2( vDirExpr)+(1-chiConv) ) );
        auto residualStabForm2 = conductivity["AIR"]*laplaciant( u )+ (gradt(u)*vDirExpr);
        form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Aq[2])
            += integrate(_range=markedelements(mesh,"AIR"),
                         _expr= coeffStab*(grad(v)*vDirExpr)*residualStabForm2 );
#endif

        // time derivative
        if ( super_type::is_time_dependent )
            for ( std::string const& marker : std::vector<std::string>({"AIR","PCB","IC1","IC2"}) )
                form2(_test=this->Xh, _trial=this->Xh, _matrix=this->M_Mq[0])
                    += integrate ( _range=markedelements( mesh,marker ),
                                   _expr=rhoC[marker]*idt( u )*id( v ) );

        // weak Dirichlet condition (rhs)
        form1(_test=this->Xh, _vector=this->M_Fq[0][0])
             += integrate( _range=markedfaces( mesh,std::list<std::string>({"Gamma_4_AIR1","Gamma_4_AIR4"} ) ),
                           _expr=300*conductivity["AIR"]*( -grad(v)*N() + penaldir*id(v)/hFace() ) );
        form1(_test=this->Xh, _vector=this->M_Fq[0][0])
            += integrate( _range=markedfaces( mesh, "Gamma_4_PCB" ),
                          _expr=300*conductivity["PCB"]*( -grad(v)*N() + penaldir*id(v)/hFace() ) );

        // Robin condition in lhs (second parameter coeff h)
        form1(_test=this->Xh, _vector=this->M_Fq[0][1])
            += integrate( _range= markedfaces( mesh,"Gamma_2" ),
                          _expr= 300*id(v) );

        // Joule effect
        form1(_test=this->Xh, _vector=this->M_Fq[0][2])
            += integrate( _range=markedelements( mesh,std::list<std::string>({"IC1"} ) ),
                          _expr=/*Qsource**/id(v) );
        form1(_test=this->Xh, _vector=this->M_Fq[0][3])
            += integrate( _range=markedelements( mesh,std::list<std::string>({"IC2"} ) ),
                          _expr=/*Qsource**/id(v) );

        // output
        form1(_test=this->Xh, _vector=this->M_Fq[1][0])
            += integrate( _range=markedelements( mesh, "IC2" ), _expr= id( v ) );
    }


#if 1
    //for scalarProduct
    auto energy = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
    energy->addMatrix(1.,this->M_Aq[0] );
    energy->addMatrix(1.,this->M_Aq[1] );
    energy->addMatrix(1.,this->M_Aq[2] );
    if ( super_type::is_time_dependent )
        if ( !boption(_name="crb.is-model-executed-in-steady-mode") )
            energy->addMatrix( M_bdf->polyDerivCoefficient( 0 ) /* 1.*/,this->M_Mq[0] );
#if 1
    energy->addMatrix(0.5,this->M_Aq[3] );
    auto uconv = this->Xh->element( vinconv*expr("-(x-0.008)*(x-0.054)*(x>0.008):x" ) );
    /*auto e = exporter( _mesh=mesh,_name="convExporter" );
     e->add( "uconv", uconv );
     e->save();*/

    auto uconvexpr = idv(uconv);
    //auto uconvexpr = vinconv*expr("-(x-0.008)*(x-0.054)*(x>0.008):x" );
    //energy->close();
    form2(_test=this->Xh, _trial=this->Xh, _matrix=energy )
        += integrate( _range= markedelements( mesh,"AIR" ),
                      _expr= 0.5*rhoC["AIR"]*(grad(u)*uconvexpr*oneY())*idt( v ) );
#endif
    energy->close();
    energy->addMatrix(1.,this->M_Aq[4] );
    this->addEnergyMatrix( energy );
#else

    auto energy2 = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
    energy2->addMatrix(1.,this->M_Aq[0] );
    energy2->addMatrix(1.,this->M_Aq[1] );
    energy2->addMatrix(1.,this->M_Aq[2] );
    energy2->close();

    auto energy = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
    energy2->symmetricPart( energy );
    energy->close();
    this->addEnergyMatrix( energy );
#endif


    //scalar product used for mass matrix
    if ( super_type::is_time_dependent )
    {
        auto InnerMassMatrix = backend()->newMatrix( _test=this->Xh, _trial=this->Xh );
        form2( _test=this->Xh, _trial=this->Xh, _matrix=InnerMassMatrix ) =
            integrate( _range=elements( mesh ), _expr=idt( u ) * id( v ) ) ;
        this->addMassMatrix(InnerMassMatrix);
    }
}


template<bool IsStationary>
double
OpusHeat<IsStationary>::output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve )
{
    using namespace vf;

    CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    double s=0;

    //bool dont_use_operators_free = boption(_name="do-not-use-operators-free") ;
    auto fqm = backend()->newVector( this->Xh );
    if ( output_index<2 )
    {
        for ( int q=0; q<M_Ql[ output_index ]; q++ )
        {
            s += this->M_betaFq[output_index][q]*dot( *this->M_Fq[output_index][q] , u );
        }
    }
    else
    {
        throw std::logic_error( "[OpusHeat::output] error with output_index : only 0 or 1 " );
    }

    return s ;
}

}


