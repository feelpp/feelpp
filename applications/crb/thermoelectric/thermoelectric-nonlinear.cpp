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

#include <feel/feelcrb/crbplugin.hpp>
#include "thermoelectric-nonlinear.hpp"

namespace Feel {

ThermoElectric::ThermoElectric()
    : super_type( "thermoelectric-nonlinear" ),
      M_propertyPath( Environment::expand( soption("thermoelectric.filename")) ),
      M_penalDir( doption("thermoelectric.penal-dir") ),
      M_picardTol( doption("thermoelectric.picard.tol") ),
      M_picardMaxit( ioption("thermoelectric.picard.maxit") ),
      M_trainsetEimSize( ioption("thermoelectric.trainset-eim-size") ),
      M_exportFE( boption("thermoelectric.export-FE") )
{}

ThermoElectric::ThermoElectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric-nonlinear" ),
      M_propertyPath( Environment::expand( soption("thermoelectric.filename")) ),
      M_penalDir( doption("thermoelectric.penal-dir") ),
      M_picardTol( doption("thermoelectric.picard.tol") ),
      M_picardMaxit( ioption("thermoelectric.picard.maxit") ),
      M_trainsetEimSize( ioption("thermoelectric.trainset-eim-size") ),
      M_exportFE( boption("thermoelectric.export-FE") )
{
    this->M_mesh = mesh;
}

/****************** Size of the decomposition *********************/
int ThermoElectric::Qa() const
{
    auto bc = M_modelProps->boundaryConditions();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    return 2*M_materials.size() + dirSize + robSize;
}

int ThermoElectric::mMaxA( int q ) const
{
    auto bc = M_modelProps->boundaryConditions();
    int matSize = M_materials.size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    std::vector<int> indexOfDirMat;
    for( auto const& exAtM : bc["potential"]["Dirichlet"])
        indexOfDirMat.push_back(indexOfMat(exAtM.material()));
    if( q < matSize )
        return this->scalarContinuousEim()[2*q]->mMax();
    else if( q < 2*matSize )
        return this->scalarContinuousEim()[2*(q-matSize)+1]->mMax();
    else if( q < 2*matSize + dirSize )
        return this->scalarContinuousEim()[2*indexOfDirMat[q-2*matSize]]->mMax();
    else if( q < 2*matSize + dirSize + robSize )
        return 1;
    else
        return 0;
}

int ThermoElectric::Nl() const
{
    return 3;
}

int ThermoElectric::Ql( int l) const
{
    auto bc = M_modelProps->boundaryConditions();
    int sourceSize = M_materials.size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    switch(l) {
    case 0:
        return sourceSize + dirSize + robSize;
    case 1:
        return 1;
    case 2:
        return 1;
    default:
        return 0;
    }
}

int ThermoElectric::mMaxL( int l, int q ) const
{
    switch( l )
    {
    case 0:
        return mMaxCompliant(q);
    case 1:
        return mMaxAverageTemp(q);
    case 2:
        return mMaxIntensity(q);
    default:
        return 0;
    }
}

int ThermoElectric::mMaxCompliant(int q ) const
{
    auto bc = M_modelProps->boundaryConditions();
    int sourceSize = M_materials.size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    std::vector<int> indexOfDirMat;
    for( auto const& exAtM : bc["potential"]["Dirichlet"])
        indexOfDirMat.push_back(indexOfMat(exAtM.material()));

    if( q < sourceSize )
        return this->scalarDiscontinuousEim()[q]->mMax();
    else if( q < sourceSize + dirSize )
        return this->scalarContinuousEim()[2*indexOfDirMat[q-sourceSize]]->mMax();
    else if( q < sourceSize + dirSize + robSize )
        return 1;
    else
        return 0;
}

int ThermoElectric::mMaxIntensity(int q ) const
{
    return 1;
}

int ThermoElectric::mMaxAverageTemp(int q ) const
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

/*************************** Initialization of the model **********************/
void ThermoElectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    M_modelProps = boost::make_shared<prop_type>(M_propertyPath);
    this->addModelFile("property-file", M_propertyPath);

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
    std::vector<std::string> range;
    M_materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    for( auto const& mat : M_materials )
        range.push_back(mat.first);

    auto domain = markedelements(M_mesh, range);
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=domain ) );

    Feel::cout << "Potential nDof  : " << Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << Xh->template functionSpace<1>()->nDof() << std::endl;

    if( !M_VT )
        M_VT = element_ptrtype( new element_type( Xh ) );
    M_V = M_VT->template elementPtr<0>();
    M_T = M_VT->template elementPtr<1>();

    M_Jh = J_space_type::New( _mesh=M_mesh, _range=domain );
    M_Th = q_sigma_space_type::New( _mesh=M_mesh, _range=domain );

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
    for( auto const& mat : M_materials )
    {
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.second.getString("sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.second.getString("alphaKey")));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto k = sigma*L*_e1;

        auto eim_sigma = eim( _model=boost::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                              _element=M_VT->template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma,
                              _space=M_Th,
                              _name=(boost::format("eim_sigma_%1%") % mat.first ).str(),
                              _sampling=Pset );
        this->addEim( eim_sigma );

        auto eim_k = eim( _model=boost::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                          _element=M_VT->template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=M_Th,
                          _name=(boost::format("eim_k_%1%") % mat.first ).str(),
                          _sampling=Pset );
        this->addEim( eim_k );

        auto eim_joule = eim( _model=boost::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                              _element=M_VT->template element<1>(),
                              _element2=M_VT->template element<0>(),
                              _parameter=M_mu,
                              _expr=sigma*_e2v*trans(_e2v),
                              _space=M_Jh,
                              _name=(boost::format("eim_joule_%1%") % mat.first ).str(),
                              _sampling=Pset );
        this->addEimDiscontinuous( eim_joule );
    }

    this->resizeQm();
    this->assemble();
}

void ThermoElectric::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    Feel::cout << "setupSpecificityModel" << std::endl;

    std::string propertyPath;
    if( this->hasModelFile("property-file") )
        propertyPath = this->additionalModelFiles().find("property-file")->second;
    else
        Feel::cerr << "Warning!! the database does not contain the property file! Expect bugs!"
                   << std::endl;
    M_modelProps = boost::make_shared<prop_type>(propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu = parameterspace_type::New( nbCrbParameters, Environment::worldComm() );

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

    this->resizeQm(false);
}

/******************************* Decomposition ***********************/
void ThermoElectric::assemble()
{
    auto VT = Xh->element();
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto phi = Xh->element();
    auto phiV = phi.template element<0>();
    auto phiT = phi.template element<1>();

    auto bc = M_modelProps->boundaryConditions();

    /************** Left hand side **************/
    int idx = 0;
    int idMat = 0;
    // electro
    for( auto const& mat : M_materials )
    {
        auto eim_sigma = this->scalarContinuousEim()[2*idMat];
        for(int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto a0 = form2(_test=Xh, _trial=Xh);
            a0 = integrate( markedelements(M_mesh, mat.first),
                            idv(eim_sigma->q(m))*inner(gradt(V),grad(phiV)) );
            M_Aqm[idx][m] = a0.matrixPtr();
        }
        idx++;
        idMat++;
    }
    idMat = 0;
    // thermo
    for( auto const& mat : M_materials )
    {
        auto eim_k = this->scalarContinuousEim()[2*idMat + 1];
        for(int m = 0; m < eim_k->mMax(); ++m )
        {
            auto a1 = form2(_test=Xh, _trial=Xh);
            a1 = integrate( markedelements(M_mesh, mat.first),
                            idv(eim_k->q(m))*inner(gradt(T),grad(phiT)) );
            M_Aqm[idx][m] = a1.matrixPtr();
        }
        idx++;
        idMat++;
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto eim_sigma = this->scalarContinuousEim()[2*indexOfMat(exAtM.material())];
        for(int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto aVD = form2(_test=Xh, _trial=Xh);
            aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                              idv(eim_sigma->q(m))*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                                    -inner(gradt(V)*N(),id(phiV))
                                                    -inner(grad(phiV)*N(),idt(V)) ) );
            M_Aqm[idx][m] = aVD.matrixPtr();
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
    for( auto const& mat : M_materials )
    {
        auto eimJoule = this->scalarDiscontinuousEim()[idx];
        for( int m = 0; m < eimJoule->mMax(); ++m )
        {
            auto f0 = form1(_test=Xh);
            f0 = integrate(elements(M_mesh),
                           inner(id(phiT), idv(eimJoule->q(m))) );
            M_Fqm[output][idx][m] = f0.vectorPtr();
        }
        idx++;
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto eim_sigma = this->scalarContinuousEim()[2*indexOfMat(exAtM.material())];
        for(int m = 0; m < eim_sigma->mMax(); ++m )
        {
            auto fVD = form1(_test=Xh);
            fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                             idv(eim_sigma->q(m))*(M_penalDir/hFace()*id(phiV) -  grad(phiV)*N()) );
            M_Fqm[output][idx][m] = fVD.vectorPtr();
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

    // Average temperature output
    auto fAvgT = form1(_test=Xh);
    for( auto const& mat : M_materials )
    {
        double area = integrate( markedelements(M_mesh, mat.first), cst(1.0) ).evaluate()(0,0) ;
        fAvgT += integrate( markedelements(M_mesh, mat.first ),
                            id(phiT)/cst(area) );
    }
    M_Fqm[output++][0][0] = fAvgT.vectorPtr();

    // Intensity output
    auto fI = form1(_test=Xh);
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            fI += integrate( markedfaces(M_mesh, exAtM.marker() ),
                             -grad(phiV)*N() );
    }
    M_Fqm[output++][0][0] = fI.vectorPtr();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(M_mesh),
                   inner(idt(V),id(phiV)) + inner(gradt(V),grad(phiV))
                   + inner(idt(T),id(phiT)) + inner(grad(T), gradt(phiT)) );
    M_M = m.matrixPtr();
    //M_energy_matrix = m.matrixPtr();

    M_InitialGuess[0][0] = Xh->elementPtr();
    M_InitialGuess[1][0] = Xh->elementPtr();
    auto TInit = M_InitialGuess[1][0]->template element<1>();
    TInit.on(elements(M_mesh), cst(1.));
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
    for( auto const& mat : M_materials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(V),grad(phiV)) );
        aqm[idx++][0] = a0.matrixPtr();
    }
    // thermo
    for( auto const& mat : M_materials )
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
    M_betaInitialGuess[1][0] = mu.parameterNamed("Tw");
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
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    for( int i = 0; i < betaEimsJoule.size(); ++i)
        betaEimsJoule[i] = this->scalarDiscontinuousEim()[i]->beta( mu, T);
    std::vector<vectorN_type> betaEimsSigma(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[2*i]->beta( mu, T);
    std::vector<vectorN_type> betaEimsK(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[2*i+1]->beta( mu, T);

    this->fillBetaQm(mu, betaEimsJoule, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    for( int i = 0; i < betaEimsJoule.size(); ++i)
        betaEimsJoule[i] = this->scalarDiscontinuousEim()[i]->beta( mu, urb);
    std::vector<vectorN_type> betaEimsSigma(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[2*i]->beta( mu, urb);
    std::vector<vectorN_type> betaEimsK(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[2*i+1]->beta( mu, urb);

    this->fillBetaQm(mu, betaEimsJoule, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    for( int i = 0; i < betaEimsJoule.size(); ++i)
        betaEimsJoule[i] = this->scalarDiscontinuousEim()[i]->beta( mu);
    std::vector<vectorN_type> betaEimsSigma(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsSigma.size(); ++i)
        betaEimsSigma[i] = this->scalarContinuousEim()[2*i]->beta( mu);
    std::vector<vectorN_type> betaEimsK(this->scalarContinuousEim().size()/2);
    for( int i = 0; i < betaEimsK.size(); ++i)
        betaEimsK[i] = this->scalarContinuousEim()[2*i+1]->beta( mu);

    this->fillBetaQm(mu, betaEimsJoule, betaEimsSigma, betaEimsK);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void ThermoElectric::fillBetaQm( parameter_type const& mu, std::vector<vectorN_type> betaEimsJoule, std::vector<vectorN_type> betaEimsSigma, std::vector<vectorN_type> betaEimsK )
{
    auto bc = M_modelProps->boundaryConditions();
    double sigmaMax = 0;
    for( auto const& mat : M_materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    // Left Hand Side
    int idx = 0;
    int idMat = 0;
    for( auto const& mat : M_materials )
    {
        for(int m = 0; m < betaEimsSigma[idMat].size(); ++m )
            M_betaAqm[idx][m] = betaEimsSigma[idMat][m]/sigmaMax;
        idx++;
        idMat++;
    }
    idMat = 0;
    for( auto const& mat : M_materials )
    {
        for( int m = 0; m < betaEimsK[idMat].size(); ++m )
            M_betaAqm[idx][m] = betaEimsK[idMat][m];
        idx++;
        idMat++;
    }
    idMat = 0;

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        idMat = indexOfMat(exAtM.material());
        for(int m = 0; m < betaEimsSigma[idMat].size(); ++m )
            M_betaAqm[idx][m] = betaEimsSigma[idMat][m]/sigmaMax;
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
    for( auto const& mat : M_materials )
    {
        for( int m = 0; m < betaEimsJoule[idx].size(); ++m )
            M_betaFqm[output][idx][m] = betaEimsJoule[idx][m];
        idx++;
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        idMat = indexOfMat(exAtM.material());
        for(int m = 0; m < betaEimsSigma[idMat].size(); ++m )
            M_betaFqm[output][idx][m] = betaEimsSigma[idMat][m]/sigmaMax*e.evaluate();
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

    // Average temperature Output
    M_betaFqm[output++][0][0] = 1;

    // Intensity Output
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        if( exAtM.expression() != "0" )
            M_betaFqm[output++][0][0] = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
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

    for( auto const& mat : M_materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        beta[idx++][0] = sigma/sigmaMax;
    }
    for( auto const& mat : M_materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        double L = mu.parameterNamed("L");
        double k = sigma*L*T;
        beta[idx++][0] = k;
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
        beta[idx++][0] = sigma/sigmaMax;
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
    for( auto const& mat : M_materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    /***************************** Picard *****************************/
    tic();
    double errV = 0, errT = 0, normOldV = 1, normOldT = 1;
    int i = 0;
    do
    {
        if( i != 0 )
        {
            normOldV = normL2(elements(M_mesh), idv(oldV));
            normOldT = normL2(elements(M_mesh), idv(oldT));
        }

        auto a = form2(_test=Xh, _trial=Xh);
        // V
        for( auto const& mat : M_materials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            a += integrate( markedelements(M_mesh, mat.first),
                            sigma/sigmaMax*inner(gradt(V),grad(phiV)) );
        }
        // T
        for( auto const& mat : M_materials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            auto L = mu.parameterNamed("L");
            auto k = sigma*cst(L)*T0;
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
        for( auto const& mat : M_materials )
        {
            auto sigma0 = mu.parameterNamed(mat.second.getString("sigmaKey"));
            auto alpha = mu.parameterNamed(mat.second.getString("alphaKey"));
            auto sigma = cst(sigma0)/(cst(1.)+cst(alpha)*(idv(oldT)-T0));
            f = integrate(markedelements(M_mesh, mat.first),
                          id(phiT)*sigma*inner(gradv(oldV), gradv(oldV)) );
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

        errV = normL2( elements(M_mesh), idv(V) - idv(oldV))/normOldV;
        errT = normL2( elements(M_mesh), idv(T) - idv(oldT))/normOldT;
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
        auto j = M_Jh->element();
        for( auto const& mat : M_materials )
        {
            auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
            j.on( markedelements(M_mesh, mat.first), _expr=sigma*inner(gradv(M_V),gradv(M_V)) );
        }
        e->add( "joule", j );
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
    if ( output_index == Nl() )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            for( int m = 0; m < mMaxL(output_index, q); ++m )
            {
                element_ptrtype eltF( new element_type( Xh ) );
                *eltF = *M_Fqm[output_index][q][0];
                output += M_betaFqm[output_index][q][0]*dot( *eltF, *M_VT );
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

int ThermoElectric::mMaxSigma()
{
    return 1;
}

ThermoElectric::q_sigma_element_type ThermoElectric::eimSigmaQ(int m)
{
    // auto Vh = Xh->template functionSpace<0>();
    q_sigma_element_type q = M_Th->element();
    q.on( _range=elements(M_mesh), _expr=cst(1.) );
    return q;
}

ThermoElectric::vectorN_type ThermoElectric::eimSigmaBeta( parameter_type const& mu )
{
    vectorN_type beta(1);
    beta(0) = mu.parameterNamed( "sigma" );
    return beta;
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

FEELPP_CRB_PLUGIN( ThermoElectric, thermoelectric )
}
