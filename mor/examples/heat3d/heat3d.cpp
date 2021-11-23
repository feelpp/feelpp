#include <feel/feelcrb/crbplugin.hpp>

#include <heat3d.hpp>

//#include <4fastsim/crbmodels/common/sourcecrbmodel.hpp>

namespace Feel
{

po::options_description
makeHeat3dOptions()
{
    po::options_description heat3doptions( "Heat3d options" );
    heat3doptions.add_options()
        ( "gamma", po::value<double>()->default_value( 10 ), "penalisation term" )
        ;
    return heat3doptions;
}
AboutData
makeHeat3dAbout( std::string const& str )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "3D Heat Application",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2016 Feel++ Consortium" );
    return about;
}


Heat3d::Heat3d()
    :
    super_type( "Heat3d" )
{}


void
Heat3d::updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const
{

}

void
Heat3d::initBetaQ()
{
    this->M_betaAq.resize( 4 );
    this->M_betaFq.resize( 2 );
    this->M_betaFq[0].resize( 4 );
    this->M_betaFq[1].resize( 1 );
}

Heat3d::super_type::betaq_type
Heat3d::computeBetaQ( parameter_type const& mu )
{
    //std::cout << "computeBetaQ start \n" << mu << std::endl;
    if ( this->M_betaAq.empty() )
        this->initBetaQ();

    for (int k = 0;k<4;++k)
        this->M_betaAq[k] = mu(k);
    for (int k = 0;k<4;++k)
        this->M_betaFq[0][k] = mu( 4+k );
    this->M_betaFq[1][0] = 1.;
    //std::cout << "computeBetaQ finish \n";
    return boost::make_tuple( this->M_betaAq, this->M_betaFq );
}


void
Heat3d::initModel()
{

    CHECK( is_linear && !is_time_dependent ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&NonLinear) == NonLinear) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";
    LOG_IF( WARNING, ((Options&TimeDependent) == TimeDependent) ) << "Invalid model is_linear:" << is_linear << " is_time_dependent:" << is_time_dependent << "\n";

    Dmu->setDimension( 8 );
    auto mu_min = Dmu->element();
    mu_min << 100,100,100,100,100*293,100*293,100*293,-1e3/*0.1*/;// fail if last param id 0 or negative
    Dmu->setMin( mu_min );
    auto mu_max = Dmu->element();
    mu_max << 200,200,200,200,200*310,200*310,200*310,1.e3;
    Dmu->setMax( mu_max );

    //return;

    auto mesh = loadMesh( _mesh=new Heat3dConfig::mesh_type,
                          _update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES),
                          _savehdf5=0 );

    this->setFunctionSpaces( Heat3dConfig::space_type::New( mesh ) );
    this->setSymbolicExpressionBuildDir("$repository/crb/heat3d/symbolicexpr/");

    if( Environment::worldComm().isMasterRank() )
    {
        std::cout << "Number of dof " << this->Xh->nDof() << "\n";
    }


    auto u = Xh->element();
    auto v = Xh->element();
    //auto mesh = Xh->mesh();

    double penaldir= doption(_name="gamma");
    //std::list<std::string> markersTimposed = {"base1","base2","base3"};
    std::vector<std::pair<std::string,std::string> > markersBases = { { "mat-base1", "base1" },{ "mat-base2", "base2" },{ "mat-base3", "base3" } };

    for (int k=0;k<3;++k )
    {
        std::string const& markerVol = markersBases[k].first;
        std::string const& markerSurf = markersBases[k].second;
        CHECK( mesh->hasMarker( markerVol ) ) << "mesh does not have the volume marker : " <<  markerVol;
        CHECK( mesh->hasMarker( markerSurf ) ) << "mesh does not have the surface marker : " <<  markerSurf;
    }
    CHECK( mesh->hasMarker( "cylinder" ) ) << "mesh does not have the surface marker : cylinder";


    double muRef = 100;

    auto energy = backend()->newMatrix(_test=Xh,_trial=Xh);

    auto a0 = form2( _trial=Xh, _test=Xh);
    a0 = integrate( _range=elements(mesh), _expr=gradt( u )*trans( grad( v ) )  );
    a0.matrixPtr()->close();

    this->addLhs( { a0 , "mu0" } );

    energy->addMatrix(muRef,a0.matrixPtr() );

    for (int k=0;k<3;++k )
    {
        std::string const& markerVol = markersBases[k].first;
        std::string const& markerSurf = markersBases[k].second;
        auto a1 = form2( _trial=Xh, _test=Xh);
        a1 = integrate( _range=markedelements(mesh,markerVol), _expr=gradt( u )*trans( grad( v ) )  );
        a1+= integrate( _range=markedfaces( mesh, markerSurf ), _expr= -(gradt(u)*id(v)+grad(v)*idt(u))*N() + penaldir*idt(u)*id(v)/hFace() );
        a1.matrixPtr()->close();

        this->addLhs( { a1 , (boost::format("mu%1%")%(k+1)).str() } );

        energy->addMatrix(muRef,a1.matrixPtr() );
    }



    int cptMark=4;
    for (int k=0;k<3;++k )
    {
        std::string const& markerSurf = markersBases[k].second;
        auto f0 = form1( _test=Xh );
        f0 = integrate( _range=markedfaces( mesh,markerSurf ), _expr=-grad(v)*N() + penaldir*id(v)/hFace() );
        f0.vectorPtr()->close();

        this->addRhs( { f0, (boost::format("mu%1%")%cptMark++).str() } );
    }

    auto f1 = form1( _test=Xh );
    f1 = integrate( _range=markedfaces( mesh,"cylinder" ), _expr=/*-*/id(v) );
    f1.vectorPtr()->close();

    this->addRhs( { f1, (boost::format("mu%1%")%cptMark++).str() } );



    /// [energy]
    //a0.matrixPtr()->symmetricPart( energy );
    //energy->addMatrix(1.,a0.matrixPtr() );
    //energy->addMatrix(1./mymu,a0.matrixPtr() );

    energy->close();
    this->addEnergyMatrix( energy );

    /// [output]
    auto out1 = form1( _test=Xh );
    double meas = integrate(_range=elements(mesh),_expr=cst(1.)).evaluate()(0,0);
    out1 = integrate( _range=elements( mesh ), _expr=id( u )/cst(meas)) ;

    this->addOutput( { out1, "1" } );

}

double
Heat3d::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve )
{
    //CHECK( ! need_to_solve ) << "The model need to have the solution to compute the output\n";

    auto mesh = Xh->mesh();
    double output=0;
    // right hand side (compliant)
    if ( output_index == 0 )
    {
        //output  = integrate( markedfaces( mesh,"BR" ), -mu(0)*expr(soption("functions.f"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0)
        //    + integrate( markedfaces( mesh,"BL" ), -mu(1)*expr(soption("functions.g"))*(gradv(u)*N()+doption("gamma")*idv(u)) ).evaluate()(0,0);
    }
    else if ( output_index == 1 )
    {
        output = mean(elements(mesh),idv(u))(0,0);
        std::cout << " Heat3d::output " << output << "\n";
    }
    // else if ( output_index == 2 )
    // {
    //     output = mean(elements(mesh),idv(u)).evaluat()(0,0);
    // }
    else
        throw std::logic_error( "[Heat3d::output] error with output_index : only 0 or 1 " );
    return output;

}

#if 0
Heat3d::parameter_type
Heat3d::crbParametersFromUserParameters( feelpp4fastsim::UserParameters const& userParam ) const
{
    parameter_type mu = Dmu->element(false);
    mu(0)=userParam[0].valueInUnitRef();
    mu(1)=userParam[1].valueInUnitRef();
    mu(2)=userParam[2].valueInUnitRef();
    mu(3)=userParam[3].valueInUnitRef();
    mu(4)=mu(1)*userParam[4].valueInUnitRef();
    mu(5)=mu(2)*userParam[5].valueInUnitRef();
    mu(6)=mu(3)*userParam[6].valueInUnitRef();
    mu(7)=userParam[7].valueInUnitRef();
    return mu;
}

void
Heat3d::updateFieldsExported( SourceCrbInterface * pvsource, element_type & feField, vtkSmartPointer<vtkUnstructuredGrid> vtkOutput )
{
    auto pvsourcecast = dynamic_cast< SourceCrbModel<Heat3d>* >( pvsource );
    if ( pvsourcecast )
    {
        pvsourcecast->addFieldExported( "Temperature", feField, vtkOutput );
    }
}
#endif

FEELPP_CRB_PLUGIN( Heat3d, heat3d )
/*
include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(heat3d)
{
    class_<Heat3d>("Heat3D")
        .def("parameterSpace", &Heat3d::)
        .def("set", &World::set)
        ;
 }*/
} // namespace Feel
