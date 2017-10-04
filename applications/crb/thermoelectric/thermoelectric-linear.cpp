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
    : super_type( "thermoelectric-linear" )
{}

ThermoElectric::ThermoElectric( mesh_ptrtype mesh )
    : super_type( "thermoelectric-linear" )
{
    this->M_mesh = mesh;
}

int ThermoElectric::Qa()
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    return 2*materials.size() + dirSize + robSize;
}

int ThermoElectric::mQA( int q )
{
    return 1;
}

int ThermoElectric::Nl()
{
    return 3;
}

int ThermoElectric::Ql( int l)
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    int sourceSize = materials.size();
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

int ThermoElectric::mLQF( int l, int q )
{
    switch( l )
    {
    case 0:
        return mCompliantQ(q);
    case 1:
        return mAverageTempQ(q);
    case 2:
        return mIntensityQ(q);
    default:
        return 0;
    }
}

int ThermoElectric::mCompliantQ(int q )
{
    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    int sourceSize = materials.size();
    int dirSize = bc["potential"]["Dirichlet"].size();
    int robSize = bc["temperature"]["Robin"].size();
    if( q < sourceSize )
        return this->scalarDiscontinuousEim()[q]->mMax();
    else if( q < sourceSize + dirSize )
        return 1;
    else if( q < sourceSize + dirSize + robSize )
        return 1;
    else
        return 0;
}

int ThermoElectric::mIntensityQ(int q )
{
    return 1;
}

int ThermoElectric::mAverageTempQ(int q )
{
    return 1;
}

void ThermoElectric::resizeQm()
{
    M_Aqm.resize( Qa());
    M_betaAqm.resize( Qa() );
    for( int q = 0; q < Qa(); ++q )
    {
        M_Aqm[q].resize(mQA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mQA(q));
    }

    M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            M_Fqm[l][q].resize(mLQF(l, q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mLQF(l, q) );
        }
    }

    M_InitialGuess.resize(1);
    M_InitialGuess[0].resize(1);
    M_InitialGuess[0][0] = Xh->elementPtr();
}

void ThermoElectric::initModel()
{
    Feel::cout << "initModel" << std::endl;
    std::string propertyPath = Environment::expand( soption("thermoelectric.filename"));
    M_modelProps = boost::make_shared<prop_type>(propertyPath);
    this->addModelFile("property-file", propertyPath);

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
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    for( auto const& mat : materials )
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

    auto Pset = this->Dmu->sampling();
    int Ne = ioption(_name="thermoelectric.trainset-eim-size");
    std::vector<size_type> N(parameterSpace()->dimension(), 1);
    N[Dmu->indexOfParameterNamed("sigma")] = Ne;
    N[Dmu->indexOfParameterNamed("current")] = Ne;

    std::string supersamplingname = (boost::format("DmuEim-Ne%1%-D%2%-generated-by-master-proc")
                                     % Ne % nbCrbParameters ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
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

    for( auto const& mat : materials )
    {
        auto eim_joule = eim( _model=boost::dynamic_pointer_cast<ThermoElectric>(this->shared_from_this() ),
                              _element=*M_V,
                              _parameter=M_mu,
                              _expr=cst_ref(M_mu.parameterNamed("sigma"))*inner(gradv(*M_V)),
                              _space=M_Jh,
                              _name=(boost::format("eim_joule_%1%") % mat.first ).str(),
                              _sampling=Pset );
        this->addEimDiscontinuous( eim_joule );
    }

    this->resizeQm();
    this->decomposition();
}

void ThermoElectric::decomposition()
{
    auto VT = Xh->element();
    auto V = VT.template element<0>();
    auto T = VT.template element<1>();
    auto phi = Xh->element();
    auto phiV = phi.template element<0>();
    auto phiT = phi.template element<1>();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    auto gamma = doption("thermoelectric.gamma");

    /************** Left hand side **************/
    int idx = 0;
    // electro
    for( auto const& mat : materials )
    {
        auto a0 = form2(_test=Xh, _trial=Xh);
        a0 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(V),grad(phiV)) );
        M_Aqm[idx++][0] = a0.matrixPtr();
    }
    // thermo
    for( auto const& mat : materials )
    {
        auto a1 = form2(_test=Xh, _trial=Xh);
        a1 = integrate( markedelements(M_mesh, mat.first),
                        inner(gradt(T),grad(phiT)) );
        M_Aqm[idx++][0] = a1.matrixPtr();
    }

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto aVD = form2(_test=Xh, _trial=Xh);
        aVD += integrate( markedfaces(M_mesh, exAtM.marker() ),
                          gamma/hFace()*inner(idt(V),id(phiV))
                          -inner(gradt(V)*N(),id(phiV))
                          -inner(grad(phiV)*N(),idt(V)) );
        M_Aqm[idx++][0] = aVD.matrixPtr();
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
    for( auto const& mat : materials )
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
        auto fVD = form1(_test=Xh);
        fVD = integrate( markedfaces(M_mesh, exAtM.marker() ),
                         gamma/hFace()*id(phiV) -  grad(phiV)*N() );
        M_Fqm[output][idx++][0] = fVD.vectorPtr();
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
    for( auto const& mat : materials )
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
                   inner(gradt(V),grad(phiV)) + inner(grad(T), gradt(phiT)) );
    M_energy_matrix = m.matrixPtr();
}

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaInitialGuess( parameter_type const& mu )
{
    M_betaInitialGuess.resize( 1 );
    M_betaInitialGuess[0].resize( 1 );
    M_betaInitialGuess[0][0] = 1;
    return this->M_betaInitialGuess;
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
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    double sigmaMax = 0;
    for( auto const& mat : materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    // Left Hand Side
    int idx = 0;
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("sigmaKey"))/sigmaMax;
    for( auto const& mat : materials )
        M_betaAqm[idx++][0] = mu.parameterNamed(mat.second.getString("kKey"));

    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
        M_betaAqm[idx++][0] = sigma/sigmaMax;
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
    for( auto const& mat : materials )
    {
        auto betaEimJoule = betaEimsJoule[idx];
        for( int m = 0; m < betaEimJoule.size(); ++m )
        {
            M_betaFqm[output][idx][m] = betaEimJoule(m);
        }
        idx++;
    }
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

        M_betaFqm[output][idx++][0] = sigma/sigmaMax*e.evaluate();
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
            M_betaFqm[output++][0][0] = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
    }
}

ThermoElectric::beta_vector_type
ThermoElectric::computeBetaLinearDecompositionA( parameter_type const& mu, double time )
{
    beta_vector_type beta;
    return beta;
}

ThermoElectric::affine_decomposition_type
ThermoElectric::computeAffineDecomposition()
{
    return boost::make_tuple( this->M_Aqm, this->M_Fqm);
}

std::vector<std::vector<ThermoElectric::sparse_matrix_ptrtype> >
ThermoElectric::computeLinearDecompositionA()
{
    return this->M_linearAqm;
}

std::vector<std::vector<ThermoElectric::element_ptrtype> >
ThermoElectric::computeInitialGuessAffineDecomposition()
{
    return M_InitialGuess;
}

ThermoElectric::element_type
ThermoElectric::solve( parameter_type const& mu )
{
    Feel::cout << "solve for parameter:" << std::endl << mu << std::endl;
    auto Vh = Xh->template functionSpace<0>();
    auto Th = Xh->template functionSpace<1>();
    auto V = Vh->element();
    auto phiV = Vh->element();
    auto T = Th->element();
    auto phiT = Th->element();

    auto bc = M_modelProps->boundaryConditions();
    auto materials = M_modelProps->materials().materialWithPhysic("thermoelectric");
    auto gamma = doption("thermoelectric.gamma");
    double sigmaMax = 0;
    for( auto const& mat : materials )
    {
        double sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        if( sigma > sigmaMax )
            sigmaMax = sigma;
    }

    /***************************** Electro *****************************/
    tic();
    auto aV = form2(_test=Vh, _trial=Vh);
    // V
    for( auto const& mat : materials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        aV += integrate( markedelements(M_mesh, mat.first),
                         sigma/sigmaMax*inner(gradt(V),grad(phiV)) );
    }

    // V Dirichlet condition
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));
        aV += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         sigma/sigmaMax*(gamma/hFace()*inner(idt(V),id(phiV))
                                         -inner(gradt(V)*N(),id(phiV))
                                         -inner(grad(phiV)*N(),idt(V)) ) );
    }

    auto fV = form1(_test=Vh);
    // V Dirichlet condition
    for( auto const& exAtM : bc["potential"]["Dirichlet"] )
    {
        auto e = expr(exAtM.expression());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );
        auto sigma = mu.parameterNamed(materials[exAtM.material()].getString("sigmaKey"));

        fV += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         sigma/sigmaMax*e.evaluate()*(gamma/hFace()*id(phiV) -  grad(phiV)*N()) );
    }
    aV.solve( _solution=M_V, _rhs=fV, _name="electro" );

    auto aT = form2(_test=Th, _trial=Th);
    // T
    for( auto const& mat : materials )
    {
        auto k = mu.parameterNamed(mat.second.getString("kKey"));
        aT += integrate( markedelements(M_mesh, mat.first),
                         k*inner(gradt(T),grad(phiT)) );
    }

    // T Robin condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression1());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        aT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         e.evaluate()*inner(idt(T), id(phiT)) );
    }

    auto fT = form1(_test=Th);
    // T right hand side
    for( auto const& mat : materials )
    {
        auto sigma = mu.parameterNamed(mat.second.getString("sigmaKey"));
        fT = integrate(markedelements(M_mesh, mat.first),
                       id(phiT)*sigma*inner(gradv(M_V), gradv(M_V)) );
    }

    // T Robin condition
    for( auto const& exAtM : bc["temperature"]["Robin"] )
    {
        auto e = expr(exAtM.expression2());
        for( auto const& param : M_modelProps->parameters() )
            if( e.expression().hasSymbol(param.first) )
                e.setParameterValues( { param.first, mu.parameterNamed(param.first) } );

        fT += integrate( markedfaces(M_mesh, exAtM.marker() ),
                         e.evaluate()*id(phiT) );
    }

    aT.solve( _solution=M_T, _rhs=fT, _name="thermo" );

    auto solution = Xh->element();
    solution.template element<0>() = M_V;
    solution.template element<1>() = M_T;

    if( boption("thermoelectric.export-FE") )
    {
        auto e = exporter(M_mesh);
        auto j = vf::project( _range=elements(M_mesh), _space=M_Jh,
                              _expr=sigma*inner(gradv(M_V),gradv(M_V)) );
        e->add( "joule", j );
        e->add("sol", solution);
        e->save();
    }
    toc("mono");

    return solution;
}

double ThermoElectric::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == Nl() )
    {
        for ( int q = 0; q < Ql(output_index); q++ )
        {
            for( int m = 0; m < mLQF(output_index, q); ++m )
            {
                element_ptrtype eltF( new element_type( Xh ) );
                *eltF = *M_Fqm[output_index][q][0];
                output += M_betaFqm[output_index][q][0]*dot( *eltF, u );
                //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
            }
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

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
    auto Vh = Xh->template functionSpace<0>();
    q_sigma_element_type q = Vh->element();
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
