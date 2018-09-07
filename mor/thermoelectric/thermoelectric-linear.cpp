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
#include "thermoelectric-linear.hpp"

namespace Feel {

ThermoElectric::ThermoElectric()
    : super_type( soption("thermoelectric.basename") ),
      M_propertyPath( Environment::expand( soption("thermoelectric.filename")) ),
      M_penalDir( doption("thermoelectric.penal-dir") ),
      M_picardTol( doption("thermoelectric.picard.tol") ),
      M_picardMaxit( ioption("thermoelectric.picard.maxit") ),
      M_trainsetEimSize( ioption("thermoelectric.trainset-eim-size") ),
      M_exportFE( boption("thermoelectric.export-FE") )
{}

ThermoElectric::ThermoElectric( mesh_ptrtype mesh )
    : super_type( soption("thermoelectric.basename") ),
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
    int elecSize = M_elecMaterials.size();
    int therSize = M_therMaterials.size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();
    return elecSize + therSize + dirVSize + dirTSize + robTSize;
}

int ThermoElectric::mMaxA( int q ) const
{
    return 1;
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
    int jouleSize = M_elecMaterials.size();
    int sourceVSize = bc["potential"]["SourceTerm"].size();
    int sourceTSize = bc["temperature"]["SourceTerm"].size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int neuVSize = bc["potential"]["Neumann"].size();
    int intVSize = bc["potential"]["Intensity"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();
    if( l == 0 )
        return jouleSize + sourceVSize + sourceTSize + dirVSize + neuVSize + intVSize + dirTSize + robTSize;
    else if( l < Nl() )
    {
        auto output = std::next(outputs.begin(), l-1)->second;
        if( output.type() == "intensity" )
            return QIntensity();
        else if( output.type() == "averageTemp")
            return QAverageTemp();
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
            return mMaxIntensity(q);
        else if( output.type() == "averageTemp")
            return mMaxAverageTemp(q);
        else
            return 0;
    }
    else
        return 0;
}

int ThermoElectric::mMaxCompliant(int q ) const
{
    auto bc = M_modelProps->boundaryConditions();
    int jouleSize = M_elecMaterials.size();
    int sourceVSize = bc["potential"]["SourceTerm"].size();
    int sourceTSize = bc["temperature"]["SourceTerm"].size();
    int dirVSize = bc["potential"]["Dirichlet"].size();
    int neuVSize = bc["potential"]["Neumann"].size();
    int intVSize = bc["potential"]["Intensity"].size();
    int dirTSize = bc["temperature"]["Dirichlet"].size();
    int robTSize = bc["temperature"]["Robin"].size();
    int index = 0;
    if( q < (index += jouleSize) )
        return this->scalarDiscontinuousEim()[0]->mMax();
    else if( q < (index += sourceVSize) )
        return 1;
    else if( q < (index += sourceTSize) )
        return 1;
    else if( q < (index += dirVSize)  )
        return 1;
    else if( q < (index += neuVSize)  )
        return 1;
    else if( q < (index += intVSize)  )
        return 1;
    else if( q < (index += dirTSize)  )
        return 1;
    else if( q < (index += robTSize)  )
        return 1;
    else
        return 0;
}

int ThermoElectric::QIntensity() const
{
    return 1;
}

int ThermoElectric::QAverageTemp() const
{
    return 1;
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
        M_initialGuess.resize(QInitialGuess());
    M_betaInitialGuess.resize(QInitialGuess());
    for( int q = 0; q < QInitialGuess(); ++q )
    {
        if( resizeMatrix )
            M_initialGuess[q].resize(mMaxInitialGuess(q), Xh->elementPtr());
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
    M_propertyPath = Environment::expand( soption("thermoelectric.filename"));
    M_penalDir = doption("thermoelectric.penal-dir");
    M_picardTol = doption("thermoelectric.picard.tol");
    M_picardMaxit = ioption("thermoelectric.picard.maxit");
    M_trainsetEimSize = ioption("thermoelectric.trainset-eim-size");
    M_exportFE = boption("thermoelectric.export-FE");

    M_modelProps = std::make_shared<prop_type>(M_propertyPath);
    M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermic"}));
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_therMaterials = M_modelProps->materials().materialWithPhysic("thermic");

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
    std::vector<bool> Nd(parameterSpace()->dimension(), false);
    // we search only for the parameters involved in the electric problem
    auto bc = M_modelProps->boundaryConditions();
    for( auto const& mat : M_elecMaterials )
        Nd[Dmu->indexOfParameterNamed(mat.second.getString("sigmaKey"))] = true;
    for( auto const& pair : bc["potential"] )
    {
        for( auto const& exAtM : pair.second )
        {
            auto e = expr(exAtM.expression1());
            for( auto const& param : M_modelProps->parameters() )
                if( e.expression().hasSymbol(param.first) )
                    Nd[Dmu->indexOfParameterNamed(param.first)] = true;
        }
    }
    for( int i = 0; i < Nd.size(); ++i )
        if( Nd[i] )
            Feel::cout << "EIM in direction of " << Dmu->parameterName(i) << std::endl;

    std::string supersamplingname = (boost::format("DmuEim-Ne%1%-D%2%-generated-by-master-proc")
                                     % M_trainsetEimSize % nbCrbParameters ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->randomInDirections( M_trainsetEimSize, Nd , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    auto gradgrad = _e2v*trans(_e2v);
    auto eim_gradgrad = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                             _element=M_VT->template element<0>(),
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=M_Jh,
                             _name="eim_gradgrad",
                             _sampling=Pset
                             );
    this->addEimDiscontinuous( eim_gradgrad );
    auto eimGradgrad = this->scalarDiscontinuousEim()[0];
    Feel::cout << "Dimension of EIM = " << eimGradgrad->mMax() << std::endl;

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
#endif
    M_mu = Dmu->element();

    if( !M_VT )
        M_VT.reset( new element_type );

    auto const& ptreeEim = ptree.get_child( "eim" );
    auto const& ptreeEimg = ptreeEim.get_child( "eim_gradgrad" );
    std::string dbnameEimg = ptreeEimg.template get<std::string>( "database-filename" );

    auto gradgrad = _e2v*trans(_e2v);
    auto eim_gradgrad = eim( _model=std::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                             _element=M_VT->template element<0>(),
                             _parameter=M_mu,
                             _expr=gradgrad,
                             _space=M_Jh,
                             _name="eim_gradgrad",
                             _filename=dbnameEimg,
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

    /************** Left hand side **************/
    int idx = 0;
    // electro
    for( auto const& mat : M_elecMaterials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 += integrate( markedelements(M_mesh, mat.first),
                         inner(gradt(V),grad(phiV)) );
        M_Aqm[idx++][0] = a0.matrixPtr();
    }
    // thermo
    for( auto const& mat : M_therMaterials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 += integrate( markedelements(M_mesh, mat.first),
                         inner(gradt(T),grad(phiT)) );
        M_Aqm[idx++][0] = a1.matrixPtr();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto aVD = form2(_test=Xh, _trial=Xh);
        aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          M_penalDir/hFace()*inner(idt(V),id(phiV))
                          -inner(gradt(V)*N(),id(phiV))
                          -inner(grad(phiV)*N(),idt(V)) );
        M_Aqm[idx++][0] = aVD.matrixPtr();
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto aTD = form2(_test=Xh, _trial=Xh);
        aTD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          M_penalDir/hFace()*inner(idt(T),id(phiT))
                          -inner(gradt(T)*N(),id(phiT))
                          -inner(grad(phiT)*N(),idt(T)) );
        M_Aqm[idx++][0] = aTD.matrixPtr();
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
    auto eimGradgrad = this->scalarDiscontinuousEim()[0];
    for( auto const& mat : M_elecMaterials )
    {
        for( int m = 0; m < eimGradgrad->mMax(); ++m )
        {
            auto f0 = form1(_test=Xh);
            f0 = integrate(markedelements(M_mesh, mat.first),
                           inner(id(phiT), idv(eimGradgrad->q(m))) );
            M_Fqm[output][idx][m] = f0.vectorPtr();
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
        auto fVD = form1(_test=Xh);
        fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         M_penalDir/hFace()*id(phiV) -  grad(phiV)*N() );
        M_Fqm[output][idx++][0] = fVD.vectorPtr();
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
        auto fTD = form1(_test=Xh);
        fTD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         M_penalDir/hFace()*id(phiT) -  grad(phiT)*N() );
        M_Fqm[output][idx++][0] = fTD.vectorPtr();
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
                auto range = markedelements(M_mesh, out.range());
                double area = integrate( range, cst(1.0) ).evaluate()(0,0) ;
                fAvgT += integrate( range, id(phiT)/cst(area) );
            }
            else if( dim == 2 )
            {
                auto range = markedfaces(M_mesh, out.range());
                double area = integrate( range, cst(1.0) ).evaluate()(0,0) ;
                fAvgT += integrate( range, id(phiT)/cst(area) );
            }
            M_Fqm[output++][0][0] = fAvgT.vectorPtr();
        }
        else if( out.type() == "intensity")
        {
            auto fI = form1(_test=Xh);
            fI += integrate( markedfaces(M_mesh, out.range() ),
                             -grad(phiV)*N() );
            M_Fqm[output++][0][0] = fI.vectorPtr();
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

    M_initialGuess[0][0] = Xh->elementPtr();
    M_initialGuess[1][0] = Xh->elementPtr();
    auto TInit = M_initialGuess[1][0]->template element<1>();
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
    Feel::cout << "computeLinearDecompositionA()" << std::endl;
    return this->M_Aqm;//M_linearAqm;
}

std::vector<std::vector<ThermoElectric::element_ptrtype> >
ThermoElectric::computeInitialGuessAffineDecomposition()
{
    return M_initialGuess;
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
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    int i = 0;
    for( auto const& eim : this->scalarDiscontinuousEim() )
        betaEimsJoule[i++] = eim->beta( mu, T);

    this->fillBetaQm(mu, betaEimsJoule);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    int i = 0;
    for( auto const& eim : this->scalarDiscontinuousEim() )
        betaEimsJoule[i++] = eim->beta( mu, urb);

    this->fillBetaQm(mu, betaEimsJoule);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

ThermoElectric::beta_type
ThermoElectric::computeBetaQm( parameter_type const& mu )
{
    std::vector<vectorN_type> betaEimsJoule(this->scalarDiscontinuousEim().size());
    int i = 0;
    for( auto const& eim : this->scalarDiscontinuousEim() )
        betaEimsJoule[i++] = eim->beta( mu);

    this->fillBetaQm(mu, betaEimsJoule);
    return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
}

void ThermoElectric::fillBetaQm( parameter_type const& mu, std::vector<vectorN_type> betaEimsJoule )
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
    for( auto const& mat : M_elecMaterials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("sigmaKey"))/sigmaMax;
    for( auto const& mat : M_therMaterials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("kKey"));

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
        M_betaAqm[idx++][0] = sigma/sigmaMax;
    }
    for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
    {
        auto k = mu.parameterNamed(M_materials[exAtM.material()].getString("kKey"));
        M_betaAqm[idx++][0] = k;
    }
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
    auto betaEimJoule = betaEimsJoule[0];
    for( auto const& mat : M_elecMaterials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        for( int m = 0; m < betaEimJoule.size(); ++m )
        {
            M_betaFqm[output][idx][m] = sigma*betaEimJoule(m);
        }
        idx++;
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
        auto sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));

        M_betaFqm[output][idx++][0] = sigma/sigmaMax*e.evaluate();
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
        auto k = mu.parameterNamed(M_materials[exAtM.material()].getString("kKey"));

        M_betaFqm[output][idx++][0] = k*e.evaluate();
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
            M_betaFqm[output++][0][0] = 1;
        }
        else if( out.type() == "intensity")
        {
            auto mat = out.getString("material");
            auto sigmaKey = M_materials[mat].getString("sigmaKey");
            M_betaFqm[output++][0][0] = mu.parameterNamed(sigmaKey);
        }
    }
}

/************************** Finite Element Resolution ********************/
ThermoElectric::element_type
ThermoElectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto VT = Xh->element();
    auto phi = Xh->element();
    auto oldVT = Xh->element();
    auto V = VT.template element<0>();
    auto phiV = phi.template element<0>();
    auto oldV = oldVT.template element<0>();
    auto T = VT.template element<1>();
    auto phiT = phi.template element<1>();
    auto oldT = oldVT.template element<1>();

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
            auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
            a += integrate( markedelements(M_mesh, mat.first),
                            sigma/sigmaMax*inner(gradt(V),grad(phiV)) );
        }
        // T
        for( auto const& mat : M_therMaterials )
        {
            auto k = mu.parameterNamed(mat.second.getString("kKey"));
            a += integrate( markedelements(M_mesh, mat.first),
                            k*inner(gradt(T),grad(phiT)) );
        }
        // V Dirichlet condition
        for( auto const& exAtM : bc["potential"]["Dirichlet"] )
        {
            auto sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));
            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            sigma/sigmaMax*(M_penalDir/hFace()*inner(idt(V),id(phiV))
                                            -inner(gradt(V)*N(),id(phiV))
                                            -inner(grad(phiV)*N(),idt(V)) ) );
        }
        // T Dirichlet condition
        for( auto const& exAtM : bc["temperature"]["Dirichlet"] )
        {
            auto k = mu.parameterNamed(M_materials[exAtM.material()].getString("kKey"));
            a += integrate( markedfaces(M_mesh, exAtM.marker() ),
                            k*(M_penalDir/hFace()*inner(idt(T),id(phiT))
                               -inner(gradt(T)*N(),id(phiT))
                               -inner(grad(phiT)*N(),idt(T)) ) );
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
            auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
            f += integrate(markedelements(M_mesh, mat.first),
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
            auto sigma = mu.parameterNamed(M_materials[exAtM.material()].getString("sigmaKey"));

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
            auto k = mu.parameterNamed(M_materials[exAtM.material()].getString("kKey"));

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
        Feel::cout << "export FE solution" << std::endl;
        auto e = exporter(M_mesh);
        // auto j = M_Jh->element();
        // for( auto const& mat : M_elecMaterials )
        // {
        //     auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        //     j.on( markedelements(M_mesh, mat.first), _expr=sigma*inner(gradv(V),gradv(V)) );
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

ThermoElectric::sparse_matrix_ptrtype ThermoElectric::energyMatrix()
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

int ThermoElectric::mMaxSigma()
{
    return M_elecMaterials.size();
}

ThermoElectric::q_sigma_element_type ThermoElectric::eimSigmaQ(int m)
{
    auto Vh = Xh->template functionSpace<0>();
    q_sigma_element_type q = Vh->element();
    auto it = M_elecMaterials.begin();
    q.on( _range=markedelements(M_mesh,std::next(it, m)->first), _expr=cst(1.) );
    return q;
}

ThermoElectric::vectorN_type ThermoElectric::eimSigmaBeta( parameter_type const& mu )
{
    vectorN_type beta(M_elecMaterials.size());
    auto it = M_elecMaterials.begin();
    for( int i = 0; i < M_elecMaterials.size(); ++i )
        beta(i) = mu.parameterNamed(std::next(it, i)->second.getString("sigmaKey"));

    return beta;
}

void ThermoElectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu )
{
    auto VT = this->solve(mu);
    auto V = VT.template element<0>();
    auto Vh = j.functionSpace();
    for( auto const& mat : M_elecMaterials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        j += vf::project(Vh, markedelements(M_mesh,mat.first), cst(-1.)*sigma*trans(gradv(V)) );
    }
}

void ThermoElectric::computeTruthCurrentDensity( current_element_type& j, parameter_type const& mu, element_type& VT )
{
    auto V = VT.template element<0>();
    auto Vh = j.functionSpace();
    for( auto const& mat : M_elecMaterials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        j += vf::project(Vh, markedelements(M_mesh,mat.first), cst(-1.)*sigma*trans(gradv(V)) );
    }
}

//FEELPP_CRB_PLUGIN( ThermoElectric, thermoelectric_linear )
FEELPP_CRBSADDLEPOINT_PLUGIN( ThermoElectric, thermoelectric_linear )
}
