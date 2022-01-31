//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#include <feel/feelmor/crbplugin.hpp>
#include "thermoelectric-nonlinear.hpp"

namespace Feel {

po::options_description
ThermoElectric::makeOptions( std::string const& prefix )
{
    po::options_description options( "Thermoelectric" );
    options.add_options()
        ( prefixvm( prefix, "basename").c_str(), Feel::po::value<std::string>()->default_value("thermoelectric_linear"),
          "name of the database" )
        ( prefixvm( prefix, "filename").c_str(), Feel::po::value<std::string>()->default_value("thermoelectric.json"),
          "json file containing application parameters and boundary conditions")
        ( prefixvm( prefix, "penal-dir").c_str(), po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( prefixvm( prefix, "trainset-eim-size").c_str(), po::value<int>()->default_value(40), "size of the eim trainset" )
        ( prefixvm( prefix, "export-FE").c_str(), po::value<bool>()->default_value(true), "export FE solution" )
        ( prefixvm( prefix, "picard.maxit").c_str(), po::value<int>()->default_value(5), "maximum number of iterations for Picard" )
        ( prefixvm( prefix, "picard.tol").c_str(), po::value<double>()->default_value(1e-8), "tolerance for Picard" )
        ;
    options.add(backend_options("thermo-electro") );
    return options;
}

ThermoElectric::ThermoElectric( std::string const& prefix )
    : super_type( soption(prefixvm( prefix, "basename").c_str()) ),
      M_propertyPath( Environment::expand( soption(prefixvm( prefix, "filename").c_str())) ),
      M_penalDir( doption(prefixvm( prefix, "penal-dir").c_str()) ),
      M_picardTol( doption(prefixvm( prefix, "picard.tol").c_str()) ),
      M_picardMaxit( ioption(prefixvm( prefix, "picard.maxit").c_str()) ),
      M_trainsetEimSize( ioption(prefixvm( prefix, "trainset-eim-size").c_str()) ),
      M_exportFE( boption(prefixvm( prefix, "export-FE").c_str()) )
{}

ThermoElectric::ThermoElectric( mesh_ptrtype mesh, std::string const& prefix )
    : super_type( soption(prefixvm( prefix, "basename").c_str()) ),
      M_propertyPath( Environment::expand( soption(prefixvm( prefix, "filename").c_str())) ),
      M_penalDir( doption(prefixvm( prefix, "penal-dir").c_str()) ),
      M_picardTol( doption(prefixvm( prefix, "picard.tol").c_str()) ),
      M_picardMaxit( ioption(prefixvm( prefix, "picard.maxit").c_str()) ),
      M_trainsetEimSize( ioption(prefixvm( prefix, "trainset-eim-size").c_str()) ),
      M_exportFE( boption(prefixvm( prefix, "export-FE").c_str()) )
{
    this->M_mesh = mesh;
}

/****************** Size of the decomposition *********************/
int ThermoElectric::Qa() const
{
    auto bc = M_modelProps->boundaryConditions();
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();
    return elecSize + therSize + dirVSize + dirTSize + robTSize;
}

int ThermoElectric::mMaxA( int q ) const
{
    auto bc = M_modelProps->boundaryConditions();
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();

    // correspondence between boundary conditions and eim
    std::vector<int> indexOfDirVMat, indexOfDirTMat;
    for( auto const& exAtM : bc["potential"]["Dirichlet"])
        indexOfDirVMat.push_back(indexOfMatV(exAtM.material()));
    for( auto const& exAtM : bc["temperature"]["Dirichlet"])
        indexOfDirTMat.push_back(indexOfMatT(exAtM.material()));

    int index = 0;
    if( q < (index += elecSize + therSize) )
        return this->scalarContinuousEim()[q]->mMax();
    else if( q < (index += dirVSize) )
        // sigma eim are the first eims
        return this->scalarContinuousEim()[indexOfDirVMat[q-(index-dirVSize)]]->mMax();
    else if( q < (index += dirTSize) )
        // k eim are after sigma eim
        return this->scalarContinuousEim()[elecSize+indexOfDirTMat[q-(index-dirTSize)]]->mMax();
    else if( q < (index += robTSize) )
        return 1;
    else
        return 0;
}

int ThermoElectric::Nl() const
{
    auto outputs = M_modelProps->outputs();
    return 1 + outputs.size();
}

int ThermoElectric::Ql( int l) const
{
    auto outputs = M_modelProps->outputs();
    auto bc = M_modelProps->boundaryConditions();
    int elecSize = M_elecMaterials.size();
    int sourceVSize = bc["potential"]["SourceTerm"].size();
    int sourceTSize = bc["temperature"]["SourceTerm"].size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int neuVSize = bc["potential"]["Neumann"].size();
    int intVSize = bc["potential"]["Intensity"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();
    if( l == 0 )
        return elecSize + sourceVSize + sourceTSize + dirVSize + neuVSize + intVSize + dirTSize + robTSize;
    else if( l < Nl() )
    {
        auto output = std::next(outputs.begin(), l-1)->second;
        if( output.type() == "intensity" )
            return QIntensity( output );
        else if( output.type() == "averageTemp")
            return QAverageTemp( output );
        else
            return 0;
    }
    else
        return 0;
}

int ThermoElectric::mMaxL( int l, int q ) const
{
    auto outputs = M_modelProps->outputs();
    if( l == 0 )
        return mMaxCompliant(q);
    else if( l < Nl() )
    {
        auto output = std::next(outputs.begin(), l-1)->second;
        if( output.type() == "intensity" )
            return mMaxIntensity(q, output);
        else if( output.type() == "averageTemp")
            return mMaxAverageTemp(q, output);
        else
            return 0;
    }
    else
        return 0;
}

int ThermoElectric::mMaxCompliant(int q ) const
{
    auto bc = M_modelProps->boundaryConditions();
    int elecSize = M_elecMaterials.size();
    int sourceVSize = bc["potential"]["SourceTerm"].size();
    int sourceTSize = bc["temperature"]["SourceTerm"].size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int neuVSize = bc["potential"]["Neumann"].size();
    int intVSize = bc["potential"]["Intensity"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();

    // correspondence between boundary conditions and eim
    std::vector<int> indexOfDirVMat, indexOfDirTMat;
    for( auto const& exAtM : bc["potential"]["Dirichlet"])
        indexOfDirVMat.push_back(indexOfMatV(exAtM.material()));
    for( auto const& exAtM : bc["temperature"]["Dirichlet"])
        indexOfDirTMat.push_back(indexOfMatT(exAtM.material()));
    int index = 0;

    if( q < (index += elecSize) )
        return this->scalarDiscontinuousEim()[0]->mMax()*this->scalarContinuousEim()[q]->mMax();
    else if( q < (index += sourceVSize) )
        return 1;
    else if( q < (index += sourceTSize) )
        return 1;
    else if( q < (index += dirVSize)  )
        return this->scalarContinuousEim()[indexOfDirVMat[q-(index-dirVSize)]]->mMax();
    else if( q < (index += neuVSize)  )
        return 1;
    else if( q < (index += intVSize)  )
        return 1;
    else if( q < (index += dirTSize)  )
        return this->scalarContinuousEim()[elecSize+indexOfDirTMat[q-(index-dirTSize)]]->mMax();
    else if( q < (index += robTSize)  )
        return 1;
    else
        return 0;
}

int ThermoElectric::QIntensity( ModelOutput const& out ) const
{
    return 1;
}

int ThermoElectric::QAverageTemp( ModelOutput const& out ) const
{
    return 1;
}

int ThermoElectric::mMaxIntensity(int q, ModelOutput const& out ) const
{
    auto mat = out.getString("material");
    auto eimSigma = this->scalarContinuousEim()[indexOfMatV(mat)];
    return eimSigma->mMax();
}

int ThermoElectric::mMaxAverageTemp(int q, ModelOutput const& out ) const
{
    return 1;
}

int ThermoElectric::QInitialGuess() const
{
    return 2;
}

int ThermoElectric::mMaxInitialGuess( int q ) const
{
    return 1;
}

void ThermoElectric::resizeQm( bool resizeMatrix )
{
    if( resizeMatrix )
        M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        if( resizeMatrix )
            M_Aqm[q].resize(mMaxA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mMaxA(q));
    }

    if( resizeMatrix )
        M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        if( resizeMatrix )
            M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            if( resizeMatrix )
                M_Fqm[l][q].resize(mMaxL(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mMaxL(l, q) );
        }
    }

    if( resizeMatrix )
        M_InitialGuess.resize(QInitialGuess());
    M_betaInitialGuess.resize(QInitialGuess());
    for( int q = 0; q < QInitialGuess(); ++q )
    {
        if( resizeMatrix )
            M_InitialGuess[q].resize(mMaxInitialGuess(q), Xh->elementPtr());
        M_betaInitialGuess[q].resize(mMaxInitialGuess(q));
    }
}

/***************************** Parameters *********************/
ThermoElectric::parameter_type
ThermoElectric::paramFromProperties() const
{
    auto mu = Dmu->element();
    int i = 0;
    auto parameters = M_modelProps->parameters();
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu(i++) = parameterPair.second.value();
        }
    }
    return mu;
}

/*************************** mesh support of functionspace **********************/
ThermoElectric::functionspace_type::mesh_support_vector_type
ThermoElectric::functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
{
    std::vector<std::string> elecRange, therRange;
    for( auto const& mat : M_elecMaterials )
        elecRange.push_back(mat.first);
    for( auto const& mat : M_therMaterials )
        therRange.push_back(mat.first);

    auto elecDomain = markedelements(mesh, elecRange);
    auto therDomain = markedelements(mesh, therRange);
    auto suppElec = std::make_shared<MeshSupport<mesh_type>>(mesh,elecDomain);
    auto suppTher = std::make_shared<MeshSupport<mesh_type>>(mesh,therDomain);
    return fusion::make_vector(suppElec,suppTher);
}

/*************************** Initialization of the model **********************/
void ThermoElectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    this->addModelFile("property-file", M_propertyPath);

    M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermic"}));
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_therMaterials = M_modelProps->materials().materialWithPhysic("thermic");

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu->setDimension(nbCrbParameters);
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    std::vector<std::string> elecRange, therRange;
    for( auto const& mat : M_elecMaterials )
        elecRange.push_back(mat.first);
    for( auto const& mat : M_therMaterials )
        therRange.push_back(mat.first);

    auto elecDomain = markedelements(M_mesh, elecRange);
    auto therDomain = markedelements(M_mesh, therRange);
    auto suppElec = std::make_shared<MeshSupport<mesh_type>>(M_mesh,elecDomain);
    auto suppTher = std::make_shared<MeshSupport<mesh_type>>(M_mesh,therDomain);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=fusion::make_vector(suppElec,suppTher) ) );

    Feel::cout << "Potential nDof  : " << Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << Xh->template functionSpace<1>()->nDof() << std::endl;

    if( !M_VT )
        M_VT = element_ptrtype( new element_type( Xh ) );
    M_V = M_VT->template elementPtr<0>();
    M_T = M_VT->template elementPtr<1>();

    M_Jh = J_space_type::New( _mesh=M_mesh, _range=elecDomain );

    auto Pset = this->Dmu->sampling();

    std::string supersamplingname = (boost::format("DmuEim-Ne%1%-D%2%-generated-by-master-proc")
                                     % M_trainsetEimSize % nbCrbParameters ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->randomize( M_trainsetEimSize , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto T0 = cst(293.0);
    for( auto const& mat : M_elecMaterials )
    {
        auto name = (boost::format("eim_sigma_%1%") % mat.first ).str();
        auto Th = q_sigma_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.first) );
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.second.getString("sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.second.getString("alphaKey")));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));

        auto eim_sigma = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                              _element=M_VT->template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma,
                              _space=Th,
                              _name=name,
                              _sampling=Pset );
        this->addEim( eim_sigma );
        Feel::cout << tc::green << name << " dimension: " << eim_sigma->mMax() << tc::reset << std::endl;
    }
    for( auto const& mat : M_therMaterials )
    {
        auto name = (boost::format("eim_k_%1%") % mat.first ).str();
        auto Th = q_sigma_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.first) );
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.second.getString("sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.second.getString("alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto k = sigma*L*_e1;

        auto eim_k = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                          _element=M_VT->template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=Th,
                          _name=name,
                          _sampling=Pset );
        this->addEim( eim_k );
        Feel::cout << tc::green << name << " dimension: " << eim_k->mMax() << tc::reset << std::endl;
    }
    auto gradgrad = _e2v*trans(_e2v);
    auto eim_gradgrad = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                             _element=M_VT->template element<0>(),
                             //_element2=M_VT->template element<0>(),
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=M_Jh,
                             _name="eim_gradgrad",
                             _sampling=Pset
                             );
    this->addEimDiscontinuous( eim_gradgrad );
    Feel::cout << tc::green << "eim_gradgrad dimension: " << eim_gradgrad->mMax() << tc::reset << std::endl;

    this->assemble();
}

void ThermoElectric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "setupSpecificityModel" << std::endl;

    if( this->hasModelFile("property-file") )
        M_propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;
    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermic"}));
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_therMaterials = M_modelProps->materials().materialWithPhysic("thermic");
#if 0
    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu = parameterspace_type::New( nbCrbParameters, Environment::worldCommPtr() );

    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& parameterPair : parameters )
    {
        if( parameterPair.second.hasMinMax() )
        {
            mu_min(i) = parameterPair.second.min();
            mu_max(i) = parameterPair.second.max();
            Dmu->setParameterName(i++, parameterPair.first );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
#endif
    M_mu = Dmu->element();

    if( !M_VT )
        M_VT.reset( new element_type );

    auto const& ptreeEim = ptree.get_child( "eim" );
    auto T0 = cst(293.0);
    for( auto const& mat : M_elecMaterials )
    {
        q_sigma_space_ptrtype Th;
        auto const& ptreeEimSigma = ptreeEim.get_child( (boost::format("eim_sigma_%1%") % mat.first ).str() );
        std::string dbnameEimSigma = ptreeEimSigma.template get<std::string>( "database-filename" );
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.second.getString("sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.second.getString("alphaKey")));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto eim_sigma = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                              _element=M_VT->template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma,
                              _space=Th,
                              _name=(boost::format("eim_sigma_%1%") % mat.first ).str(),
                              _filename=dbnameEimSigma,
                              _directory=dbDir );
        this->addEim( eim_sigma );
    }
    for( auto const& mat : M_therMaterials )
    {
        q_sigma_space_ptrtype Th;
        auto const& ptreeEimK = ptreeEim.get_child( (boost::format("eim_k_%1%") % mat.first ).str() );
        std::string dbnameEimK = ptreeEimK.template get<std::string>( "database-filename" );
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.second.getString("sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.second.getString("alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto k = sigma*L*_e1;
        auto eim_k = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                          _element=M_VT->template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=Th,
                          _name=(boost::format("eim_k_%1%") % mat.first ).str(),
                          _filename=dbnameEimK,
                          _directory=dbDir );
        this->addEim( eim_k );
    }
    auto const& ptreeEimGrad = ptreeEim.get_child( "eim_gradgrad" );
    std::string dbnameEimGrad = ptreeEimGrad.template get<std::string>( "database-filename" );

    auto gradgrad = _e2v*trans(_e2v);
    auto eim_gradgrad = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                             _element=M_VT->template element<0>(),
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=M_Jh,
                             _name="eim_gradgrad",
                             _filename=dbnameEimGrad,
                             _directory=dbDir );
    this->addEimDiscontinuous( eim_gradgrad );

    this->resizeQm(false);
}

/******************************* Decomposition ***********************/
void ThermoElectric::assemble()
{
    this->resizeQm();

    auto VT = Xh->element();
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto phi = Xh->element();
    auto phiV = phi.template element<0>();
    auto phiT = phi.template element<1>();

    auto bc = M_modelProps->boundaryConditions();
    int elecSize = M_elecMaterials.size();

    /************** Left hand side **************/
    int idx = 0;
    int idMat = 0;
    // electro
    for( auto const& mat : M_elecMaterials )
    {
        auto eimSigma = this->scalarContinuousEim()[idMat];
        for(int m = 0; m < eimSigma->mMax(); ++m )
        {
            auto a0 = form2(_test=Xh, _trial=Xh);
            a0 = integrate( markedelements(M_mesh, mat.first),
                            idv(eimSigma->q(m))*inner(gradt(V),grad(phiV)) );
            M_Aqm[idx][m] = a0.matrixPtr();
        }
        idx++;
        idMat++;
    }
    idMat = 0;
    // thermo
    for( auto const& mat : M_therMaterials )
    {
        auto eimK = this->scalarContinuousEim()[idMat + elecSize];
        for(int m = 0; m < eimK->mMax(); ++m )
        {
            auto a1 = form2(_test=Xh, _trial=Xh);
            a1 = integrate( markedelements(M_mesh, mat.first),
                            idv(eimK->q(m))*inner(gradt(T),grad(phiT)) );
            M_Aqm[idx][m] = a1.matrixPtr();
        }
        idx++;
        idMat++;
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto eimSigma = this->scalarContinuousEim()[indexOfMatV(exAtM.material())];
        for(int m = 0; m < eimSigma->mMax(); ++m )
        {
            auto aVD = form2(_test=Xh, _trial=Xh);
            aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                              idv(eimSigma->q(m))*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                                    -inner(gradt(V)*N(),id(phiV))
                                                    -inner(grad(phiV)*N(),idt(V)) ) );
            M_Aqm[idx][m] = aVD.matrixPtr();
        }
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto eimK = this->scalarContinuousEim()[elecSize+indexOfMatT(exAtM.material())];
        for(int m = 0; m < eimK->mMax(); ++m )
        {
            auto aTD = form2(_test=Xh, _trial=Xh);
            aTD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                              idv(eimK->q(m))*(M_penalDir/hFace()*inner(idt(T),id(phiT))
                                                -inner(gradt(T)*N(),id(phiT))
                                                -inner(grad(phiT)*N(),idt(T)) ) );
            M_Aqm[idx][m] = aTD.matrixPtr();
        }
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto aTR = form2(_test=Xh, _trial=Xh);
        aTR += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          inner(idt(T),id(phiT)) );
        M_Aqm[idx++][0] = aTR.matrixPtr();
    }

    /************** Right hand side **************/
    int output = 0;
    idx = 0;
    auto eimGrad = this->scalarDiscontinuousEim()[0];
    for( auto const& mat : M_elecMaterials )
    {
        auto eimSigma = this->scalarContinuousEim()[idx];
        int mSigma = eimSigma->mMax();
        for( int k = 0; k < eimGrad->mMax(); ++k )
        {
            auto gradk = eimGrad->q(k);
            for( int m = 0; m < mSigma; ++m )
            {
                auto fJ = form1(_test=Xh);
                fJ = integrate(markedelements(M_mesh, mat.first),
                               idv(eimSigma->q(m))*inner(id(phiT), idv(gradk)) );
                M_Fqm[output][idx][k*mSigma+m] = fJ.vectorPtr();
            }
        }
        idx++;
    }
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression2());
        auto fVST = form1(_test=Xh);
        fVST = integrate( markedelements(M_mesh, exAtM.marker() ),
                          e*id(phiV) );
        M_Fqm[output][idx++][0] = fVST.vectorPtr();
    }
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression2());
        auto fTST = form1(_test=Xh);
        fTST = integrate( markedelements(M_mesh, exAtM.marker() ),
                          e*id(phiT) );
        M_Fqm[output][idx++][0] = fTST.vectorPtr();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto eimSigma = this->scalarContinuousEim()[indexOfMatV(exAtM.material())];
        for(int m = 0; m < eimSigma->mMax(); ++m )
        {
            auto fVD = form1(_test=Xh);
            fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                             idv(eimSigma->q(m))*(M_penalDir/hFace()*id(phiV) -  grad(phiV)*N()) );
            M_Fqm[output][idx][m] = fVD.vectorPtr();
        }
        idx++;
    }
    for( auto const& exAtM : bc["potential"]["Neumann"] )
    {
        auto fVN = form1(_test=Xh);
        fVN = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         id(phiV) );
        M_Fqm[output][idx++][0] = fVN.vectorPtr();
    }
    for( auto const& exAtM : bc["potential"]["Intensity"] )
    {
        double area = integrate( markedfaces(M_mesh, exAtM.marker() ), cst(1.) ).evaluate()(0,0);
        auto fVI = form1(_test=Xh);
        fVI = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         id(phiV)/area );
        M_Fqm[output][idx++][0] = fVI.vectorPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto eimK = this->scalarContinuousEim()[elecSize+indexOfMatT(exAtM.material())];
        for(int m = 0; m < eimK->mMax(); ++m )
        {
            auto fTD = form1(_test=Xh);
            fTD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                             idv(eimK->q(m))*(M_penalDir/hFace()*id(phiT) -  grad(phiT)*N()) );
            M_Fqm[output][idx][m] = fTD.vectorPtr();
        }
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto fTR = form1(_test=Xh);
        fTR = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         id(phiT) );
        M_Fqm[output][idx++][0] = fTR.vectorPtr();
    }
    output++;

    auto outputs = M_modelProps->outputs();
    for( auto const& outp : outputs )
    {
        auto out = outp.second;
        if( out.type() == "averageTemp" )
        {
            auto dim = out.dim();
            auto fAvgT = form1(_test=Xh);
            if( dim == 3 )
            {
                auto range = markedelements(M_mesh, out.markers());
                double area = integrate( range, cst(1.0) ).evaluate()(0,0) ;
                fAvgT += integrate( range, id(phiT)/cst(area) );
            }
            else if( dim == 2 )
            {
                auto range = markedfaces(M_mesh, out.markers());
                double area = integrate( range, cst(1.0) ).evaluate()(0,0) ;
                fAvgT += integrate( range, id(phiT)/cst(area) );
            }
            M_Fqm[output++][0][0] = fAvgT.vectorPtr();
        }
        else if( out.type() == "intensity")
        {
            auto mat = out.getString("material");
            auto eimSigma = this->scalarContinuousEim()[indexOfMatV(mat)];
            for( int m = 0; m < eimSigma->mMax(); ++m )
            {
                auto fI = form1(_test=Xh);
                fI += integrate( markedfaces(M_mesh, out.markers() ),
                                 -idv(eimSigma->q(m))*grad(phiV)*N() );
                M_Fqm[output][0][m] = fI.vectorPtr();
            }
            output++;
        }
    }

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    for( auto const& mat : M_elecMaterials )
    {
        m += integrate( markedelements(M_mesh, mat.first),
                        inner(idt(V),id(phiV)) + inner(gradt(V),grad(phiV)) );
    }
    for( auto const& mat : M_therMaterials )
    {
        m += integrate( markedelements(M_mesh, mat.first),
                        inner(idt(T),id(phiT)) + inner(grad(T), gradt(phiT)) );
    }
    M_M = m.matrixPtr();
    M_energy_matrix = m.matrixPtr();

    M_InitialGuess[0][0] = Xh->elementPtr();
    M_InitialGuess[1][0] = Xh->elementPtr();
    auto TInit = M_InitialGuess[1][0]->template element<1>();
    for( auto const& mat : M_materials )
        TInit.on(markedelements(M_mesh,mat.first), cst(293.));
}

ThermoElectric::affine_decomposition_type
ThermoElectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<ThermoElectric::sparse_matrix_ptrtype> >
ThermoElectric::computeLinearDecompositionA()
{
    auto aqm = std::vector<std::vector<sparse_matrix_ptrtype> >();
    aqm.resize(Qa(), std::vector<sparse_matrix_ptrtype>(1));
    auto VT = Xh->element();
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto phi = Xh->element();
    auto phiV = phi.template element<0>();
    auto phiT = phi.template element<1>();

    auto bc = M_modelProps->boundaryConditions();

    /************** Left hand side **************/
    int idx = 0;
    // electro
    for( auto const& mat : M_elecMaterials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(V),grad(phiV)) );
        aqm[idx++][0] = a0.matrixPtr();
    }
    // thermo
    for( auto const& mat : M_therMaterials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(T),grad(phiT)) );
        aqm[idx++][0] = a1.matrixPtr();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto aVD = form2(_test=Xh, _trial=Xh);
        aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          M_penalDir/hFace()*inner(idt(V),id(phiV))
                          -inner(gradt(V)*N(),id(phiV))
                          -inner(grad(phiV)*N(),idt(V)) );
        aqm[idx++][0] = aVD.matrixPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto aTD = form2(_test=Xh, _trial=Xh);
        aTD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          M_penalDir/hFace()*inner(idt(T),id(phiT))
                          -inner(gradt(T)*N(),id(phiT))
                          -inner(grad(phiT)*N(),idt(T)) );
        aqm[idx++][0] = aTD.matrixPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto aTR = form2(_test=Xh, _trial=Xh);
        aTR += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          inner(idt(T),id(phiT)) );
        aqm[idx++][0] = aTR.matrixPtr();
    }

    return aqm;
}

std::vector<std::vector<ThermoElectric::element_ptrtype> >
ThermoElectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

/***************************** Beta coefficients ****************************/

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess[0][0] = 0;
    M_betaInitialGuess[1][0] = 1;
    return this->M_betaInitialGuess;
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( parameter_type const& mu, double time, bool only_terms_time_dependent )
{
    return computeBetaQm(mu);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( element_type const& T, parameter_type const& mu )
{
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta( mu, T);
    std::vector<vectorN_type> betaEimsSigma(elecSize);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[i]->beta( mu, T);
    std::vector<vectorN_type> betaEimsK(therSize);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[i+elecSize]->beta( mu, T);

    this->fillBetaQm(mu, betaEimGrad, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta( mu, urb);
    std::vector<vectorN_type> betaEimsSigma(elecSize);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[i]->beta( mu, urb);
    std::vector<vectorN_type> betaEimsK(therSize);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[i+elecSize]->beta( mu, urb);

    this->fillBetaQm(mu, betaEimGrad, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( parameter_type const& mu )
{
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    vectorN_type betaEimGrad = this->scalarDiscontinuousEim()[0]->beta( mu );
    std::vector<vectorN_type> betaEimsSigma(elecSize);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[i]->beta( mu);
    std::vector<vectorN_type> betaEimsK(therSize);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[i+elecSize]->beta( mu);

    this->fillBetaQm(mu, betaEimGrad, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void ThermoElectric::fillBetaQm( parameter_type const& mu, vectorN_type betaEimGrad, std::vector<vectorN_type> betaEimsSigma, std::vector<vectorN_type> betaEimsK )
{
    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& mat : M_elecMaterials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    // Left Hand Side
    int idx = 0;
    int idMat = 0;
    for( auto const& mat : M_elecMaterials )
    {
        int mMax = std::min<int>(betaEimsSigma[idMat].size(),  M_betaAqm[idx].size() );
        for(int m = 0; m < mMax; ++m )
            M_betaAqm[idx][m] = betaEimsSigma[idMat][m]/sigmaMax;
        idx++;
        idMat++;
    }
    idMat = 0;
    for( auto const& mat : M_therMaterials )
    {
        int mMax = std::min<int>(betaEimsK[idMat].size(), M_betaAqm[idx].size());
        for(int m = 0; m < mMax; ++m )
            M_betaAqm[idx][m] = betaEimsK[idMat][m];
        idx++;
        idMat++;
    }

    idMat = 0;
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        idMat = indexOfMatV(exAtM.material());
        int mMax = std::min<int>(betaEimsSigma[idMat].size(), M_betaAqm[idx].size());
        for(int m = 0; m < mMax; ++m )
            M_betaAqm[idx][m] = betaEimsSigma[idMat][m]/sigmaMax;
        idx++;
    }
    idMat = 0;
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        idMat = indexOfMatT(exAtM.material());
        int mMax = std::min<int>(betaEimsK[idMat].size(), M_betaAqm[idx].size());
        for(int m = 0; m < mMax; ++m )
        for(int m = 0; m < betaEimsK[idMat].size(); ++m )
            M_betaAqm[idx][m] = betaEimsK[idMat][m];
        idx++;
    }
    idMat = 0;
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaAqm[idx++][0] = e.evaluate();
    }

    // Right Hand Side
    int output = 0;
    idx = 0;
    idMat = 0;
    for( auto const& mat : M_elecMaterials )
    {
        int mGrad = betaEimGrad.size();
        int mSigma = betaEimsSigma[idMat].size() < M_betaFqm[output][idx].size()/mGrad? betaEimsSigma[idMat].size() : M_betaFqm[output][idx].size()/mGrad;
        for( int k = 0; k < mGrad; ++k )
        {
            for( int m = 0; m < mSigma; ++m )
                M_betaFqm[output][idx][k*mSigma+m] = betaEimGrad[k]*betaEimsSigma[idMat][m];
        }
        idx++;
        idMat++;
    }
    for( auto const& exAtM : bc["potential"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate();
    }
    for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        idMat = indexOfMatV(exAtM.material());
        int mMax = std::min<int>(betaEimsSigma[idMat].size(), M_betaFqm[output][idx].size());
        for(int m = 0; m < mMax; ++m )
            M_betaFqm[output][idx][m] = betaEimsSigma[idMat][m]/sigmaMax*e.evaluate();
        idx++;
    }
    for( auto const& exAtM : bc["potential"]["Neumann"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate()/sigmaMax;
    }
    for( auto const& exAtM : bc["potential"]["Intensity"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate()/sigmaMax;
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        idMat = indexOfMatT(exAtM.material());
        int mMax = std::min<int>(betaEimsK[idMat].size(), M_betaFqm[output][idx].size());
        for(int m = 0; m < mMax; ++m )
            M_betaFqm[output][idx][m] = betaEimsK[idMat][m]*e.evaluate();
        idx++;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        M_betaFqm[output][idx++][0] = e.evaluate();
    }
    output++;

    auto outputs = M_modelProps->outputs();
    for( auto const& outp : outputs )
    {
        auto out = outp.second;
        if( out.type() == "averageTemp" )
        {
            M_betaFqm[output][0][0] = 1;
        }
        else if( out.type() == "intensity")
        {
            idMat = indexOfMatV(out.getString("material"));
            int mMax = std::min<int>(betaEimsSigma[idMat].size(), M_betaFqm[output][0].size());
            for( int m = 0; m < mMax; ++m )
                M_betaFqm[output][0][m] = betaEimsSigma[idMat][m];
        }
        output++;
    }
}

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaLinearDecompositionA( parameter_type const& mu , double time )
{
    beta_vector_type beta(Qa(), std::vector<double>(1));
    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& mat : M_materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }
    double T = 293.;

    int idx = 0;

    for( auto const& mat : M_elecMaterials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        beta[idx++][0] = sigma/sigmaMax;
    }
    for( auto const& mat : M_therMaterials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        double L = mu.parameterNamed("L");
        double k = sigma*L*T;
        beta[idx++][0] = k;
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        double sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
        beta[idx++][0] = sigma/sigmaMax;
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        double sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
        double L = mu.parameterNamed("L");
        double k = sigma*L*T;
        beta[idx++][0] = k;
    }
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        beta[idx++][0] = e.evaluate();
    }

    return beta;
}

/************************** Finite Element Resolution ********************/
ThermoElectric::element_type
ThermoElectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto T0 = cst(293.0);

    auto VT = Xh->element();
    auto phi = Xh->element();
    auto oldVT = Xh->element();
    auto V = VT.template element<0>();
    auto phiV = phi.template element<0>();
    auto oldV = oldVT.template element<0>();
    auto T = VT.template element<1>();
    auto phiT = phi.template element<1>();
    auto oldT = oldVT.template element<1>( T0 ); // init T with 293

    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& mat : M_elecMaterials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    /***************************** Picard *****************************/
    tic();
    double errV = 0, errT = 0, errOldV = 1, errOldT = 1;
    int i = 0;
    do
    {
        if( i != 0 )
        {
            errOldV = 0;
            errOldT = 0;
            for( auto const& mat : M_elecMaterials )
                errOldV += normL2(markedelements(M_mesh,mat.first), idv(oldV));
            for( auto const& mat : M_therMaterials )
                errOldT += normL2(markedelements(M_mesh,mat.first), idv(oldT));
        }

        auto a = form2(_test=Xh, _trial=Xh);
        // V
        for( auto const& mat : M_elecMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            a += integrate( markedelements(M_mesh, mat.first),
                            sigma/sigmaMax*inner(gradt(V),grad(phiV)) );
        }
        // T
        for( auto const& mat : M_therMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            auto L = mu.parameterNamed("L");
            auto k = sigma*cst(L)*idv(oldT);
            a += integrate( markedelements(M_mesh, mat.first),
                            k*inner(gradt(T),grad(phiT)) );
        }
        // V Dirichlet condition
        for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        {
            auto sigma0 = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
            auto alpha = mu.parameterNamed(M_materials[exAtM.material()].getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            sigma/sigmaMax*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                            -inner(gradt(V)*N(),id(phiV))
                                            -inner(grad(phiV)*N(),idt(V)) ) );
        }
        // T Dirichlet condition
        for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
        {
            auto sigma0 = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
            auto alpha = mu.parameterNamed(M_materials[exAtM.material()].getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            auto L = mu.parameterNamed("L");
            auto k = sigma*cst(L)*idv(oldT);
            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            k*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                               -inner(gradt(V)*N(),id(phiV))
                               -inner(grad(phiV)*N(),idt(V)) ) );
        }
        // T Robin condition
        for( auto const& exAtM : bc["temperature"]["Robin"] )
        {
            auto e = expr(exAtM.expression1());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            e.evaluate()*inner(idt(T), id(phiT)) );
        }

        auto f = form1(_test=Xh);
        // T right hand side
        for( auto const& mat : M_elecMaterials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            f = integrate(markedelements(M_mesh, mat.first),
                          id(phiT)*sigma*inner(gradv(oldV), gradv(oldV)) );
        }
        // V source term
        for( auto const& exAtM : bc["potential"]["SourceTerm"] )
        {
            auto e1 = expr(exAtM.expression1());
            for( auto const& param : M_modelProps->parameters() )
                if( e1.expression().hasSymbol(param.first) )
                    e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            auto e2 = expr(exAtM.expression2());

            f += integrate( markedelements(M_mesh, exAtM.marker() ),
                            e2*e1.evaluate()*id(phiV)/sigmaMax );
        }
        // T source term
        for( auto const& exAtM : bc["temperature"]["SourceTerm"] )
        {
            auto e1 = expr(exAtM.expression1());
            for( auto const& param : M_modelProps->parameters() )
                if( e1.expression().hasSymbol(param.first) )
                    e1.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            auto e2 = expr(exAtM.expression2());

            f += integrate( markedelements(M_mesh, exAtM.marker() ),
                            e2*e1.evaluate()*id(phiT) );
        }
        // V Dirichlet condition
        for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        {
            auto e = expr(exAtM.expression());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            auto sigma0 = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
            auto alpha = mu.parameterNamed(M_materials[exAtM.material()].getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            sigma/sigmaMax*e.evaluate()*(M_penalDir/hFace()*id(phiV) -  grad(phiV)*N()) );
        }
        // V Neumann condition
        for( auto const& exAtM : bc["potential"]["Neumann"] )
        {
            auto e = expr(exAtM.expression());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            e.evaluate()*id(phiV)/sigmaMax );
        }
        // V Intensity condition
        for( auto const& exAtM : bc["potential"]["Intensity"] )
        {
            double area = integrate( markedfaces(M_mesh, exAtM.marker()), cst(1.) ).evaluate()(0,0);
            auto e = expr(exAtM.expression());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            e.evaluate()/area*id(phiV)/sigmaMax );
        }
        // T Dirichlet condition
        for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
        {
            auto e = expr(exAtM.expression());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
            auto sigma0 = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
            auto alpha = mu.parameterNamed(M_materials[exAtM.material()].getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            auto L = mu.parameterNamed("L");
            auto k = sigma*cst(L)*idv(oldT);

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            k*e.evaluate()*(M_penalDir/hFace()*id(phiT) -  grad(phiT)*N()) );
        }
        // T Robin condition
        for( auto const& exAtM : bc["temperature"]["Robin"] )
        {
            auto e = expr(exAtM.expression2());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

            f += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            e.evaluate()*id(phiT) );
        }

        a.solve( _solution=VT, _rhs=f, _name="thermo-electro" );

        errV = 0;
        errT = 0;
        for( auto const& mat : M_elecMaterials )
            errV += normL2( markedelements(M_mesh, mat.first), idv(V) - idv(oldV));
        for( auto const& mat : M_therMaterials )
            errT += normL2( markedelements(M_mesh, mat.first), idv(T) - idv(oldT));
        errV /= errOldV;
        errT /= errOldT;
        Feel::cout << "Picard[" << i << "] err(V) = " << errV << " err(T) = " << errT << std::endl;

        oldV = V;
        oldT = T;
    }
    while( ++i < M_picardMaxit && (errT > M_picardTol || errV > M_picardTol) );

    auto solution = Xh->element();
    solution.template element<0>() = V;
    solution.template element<1>() = T;

    if( M_exportFE )
    {
        auto e = exporter(M_mesh);
        // auto j = M_Jh->element();
        // for( auto const& mat : M_materials )
        // {
        //     auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        //     j.on( markedelements(M_mesh, mat.first), _expr=sigma*inner(gradv(M_V),gradv(M_V)) );
        // }
        // e->add( "joule", j );
        e->add("sol", solution);
        e->save();
    }
    toc("mono");

    return solution;
}

/**************************** Scalar product ********************/
double
ThermoElectric::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M_M->energy( x, y );
}

double
ThermoElectric::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M_M->energy( x, y );
}

ThermoElectric::sparse_matrix_ptrtype
ThermoElectric::energyMatrix()
{
    return M_M;
}


/**************************** Output *********************/
double ThermoElectric::output( int output_index, parameter_type const& mu , element_type& solution, bool need_to_solve)
{
    Feel::cout << "compute output " << output_index << " for mu=\n" << mu << std::endl;
    this->computeBetaQm( solution, mu );

    if ( need_to_solve )
    {
        *M_VT = this->solve( mu );
        //this->update( mu );
        //backend()->solve( _matrix=D,  _solution=pT, _rhs=F );
    }
    else
        *M_VT = solution;

    double output=0;
    if ( output_index < Nl() )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            for( int m = 0; m < mMaxL(output_index, q); ++m )
            {
                element_ptrtype eltF( new element_type( Xh ) );
                *eltF = *M_Fqm[output_index][q][m];
                output += M_betaFqm[output_index][q][m]*dot( *eltF, *M_VT );
                //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
            }
        }
    }
    else
        throw std::logic_error( "[Thermoelectric::output] error with output_index :\n 0: compliant, 1: temperature average, 2: intensity" );
    return output;
}

/******************************* BiotSavart API ************************/
int ThermoElectric::mMaxJoule()
{
    return this->scalarDiscontinuousEim()[0]->mMax();
}

int ThermoElectric::mMaxSigma( std::string const& mat )
{
    return this->scalarContinuousEim()[indexOfMatV(mat)]->mMax();
}

ThermoElectric::q_sigma_element_type ThermoElectric::eimSigmaQ( std::string const& mat, int m)
{
    return this->scalarContinuousEim()[indexOfMatV(mat)]->q(m);
}

ThermoElectric::vectorN_type ThermoElectric::eimSigmaBeta( std::string const& mat, parameter_type const& mu )
{
    return this->scalarContinuousEim()[indexOfMatV(mat)]->beta(mu);
}

ThermoElectric::vectorN_type ThermoElectric::eimSigmaBeta( std::string const& mat, parameter_type const& mu, element_type const& U )
{
    return this->scalarContinuousEim()[indexOfMatV(mat)]->beta(mu, U);
}

ThermoElectric::vectorN_type ThermoElectric::eimSigmaBeta( std::string const& mat, parameter_type const& mu, vectorN_type const& Urb )
{
    return this->scalarContinuousEim()[indexOfMatV(mat)]->beta(mu, Urb);
}

void ThermoElectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu )
{
    auto VT = this->solve(mu);
    auto V = VT.template element<0>();
    auto sigma = mu.parameterNamed("sigma");
    auto Vh = j.functionSpace();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
}

void ThermoElectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu, element_type& VT )
{
    auto V = VT.template element<0>();
    auto sigma = mu.parameterNamed("sigma");
    auto Vh = j.functionSpace();
    j = vf::project(Vh, elements(M_mesh), cst(-1.)*sigma*trans(gradv(V)) );
}

FEELPP_CRBSADDLEPOINT_PLUGIN( ThermoElectric, thermoelectric )
}
