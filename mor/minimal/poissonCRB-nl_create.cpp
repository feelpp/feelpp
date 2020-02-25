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
#include "poissonCRB-nl.hpp"

namespace Feel {

po::options_description
PoissonNL::makeOptions()
{
    FEELPP_EXPORT po::options_description options( "Poisson" );
    options.add_options()
        ( "poisson.filename", Feel::po::value<std::string>()->default_value("poisson.json"),
          "json file containing application parameters")
        ( "poisson.gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( "poisson.tolerance", po::value<double>()->default_value(1e-8), "tolerance for picard" )
        ( "poisson.maxit", po::value<int>()->default_value(50), "maximum number of iteration for picard" )
        ( "poisson.trainset-eim-size", po::value<int>()->default_value(10), "size of the trainset" )
        ;
    options.add(backend_options("thermo-electro") );
    options.add(backend_options("electro") );
    options.add(backend_options("thermo") );
    return options;
}

AboutData
PoissonNL::makeAbout( std::string const& str )
{
    AboutData about( str.c_str() );
    return about;
}



PoissonNL::PoissonNL() : super_type("poissonmodel-nl_crb") {}

int PoissonNL::Qa()
{
    return M_nbTherMat + M_nbElecMat + M_nbPotDir + M_nbTempRobin;
}

int PoissonNL::mMaxA(int q)
{
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    // correspondence between boundary conditions and eim
    std::vector<int> indexOfPotDirMat;
    for( auto const&[bdName, bd] : potentialDirichlet )
        indexOfPotDirMat.push_back(indexOfElecMat(bd.material()));

    if( q < M_nbTherMat )
        return this->scalarContinuousEim()[q]->mMax();
    if( q < M_nbTherMat + M_nbElecMat )
        return this->scalarContinuousEim()[q]->mMax();
    if( q < M_nbTherMat + M_nbElecMat + M_nbPotDir )
        return this->scalarContinuousEim()[indexOfPotDirMat[q-(M_nbTherMat+M_nbElecMat)]+M_nbTherMat]->mMax();
    if( q < M_nbTherMat + M_nbElecMat + M_nbPotDir + M_nbTempRobin )
        return 1;
    return 0;
}

int PoissonNL::Nl()
{
    return 1;
}

int PoissonNL::Ql( int l)
{
    switch(l) {
    case 0:
        return M_nbElecMat + M_nbPotDir + M_nbTempRobin;
    default:
        return 0;
    }
}

int PoissonNL::mMaxF(int l, int q)
{
    if( l == 0 )
    {
        auto bdConditions = M_modelProps->boundaryConditions2();
        auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
        // correspondence between boundary conditions and eim
        std::vector<int> indexOfPotDirMat;
        for( auto const&[bdName, bd] : potentialDirichlet )
            indexOfPotDirMat.push_back(indexOfElecMat(bd.material()));

        if( q < M_nbElecMat )
            return this->scalarContinuousEim()[q+M_nbTherMat]->mMax()*this->scalarDiscontinuousEim()[0]->mMax();
        if( q < M_nbElecMat + M_nbPotDir )
            return this->scalarContinuousEim()[indexOfPotDirMat[q-M_nbElecMat]+M_nbTherMat]->mMax();
        if( q < M_nbElecMat + M_nbPotDir + M_nbTempRobin )
            return 1;
        return 0;
    }
    else
        return 0;
}

void PoissonNL::resize()
{
    M_Aqm.resize(Qa());
    M_betaAqm.resize(Qa());
    for( int q = 0; q < Qa(); ++q)
    {
        M_Aqm[q].resize(mMaxA(q), backend()->newMatrix(Xh, Xh ) );
        M_betaAqm[q].resize(mMaxA(q));
    }
    M_Fqm.resize(Nl());
    M_betaFqm.resize(Nl());
    for( int l = 0; l < Nl(); ++l )
    {
        M_Fqm[l].resize(Ql(l));
        M_betaFqm[l].resize(Ql(l));
        for( int q = 0; q < Ql(l); ++q )
        {
            M_Fqm[l][q].resize(mMaxF(l,q), backend()->newVector(Xh) );
            M_betaFqm[l][q].resize(mMaxF(l,q));
        }
    }
}

int PoissonNL::indexOfElecMat(std::string const& mat ) const
{
    auto it = M_elecMaterials.find(mat);
    if( it != M_elecMaterials.end() )
        return std::distance(it, M_elecMaterials.begin());
    else
        return -1;
}

int PoissonNL::indexOfTherMat(std::string const& mat ) const
{
    auto it = M_therMaterials.find(mat);
    if( it != M_therMaterials.end() )
        return std::distance(it, M_therMaterials.begin());
    else
        return -1;
}


/*************************** mesh support of functionspace **********************/
PoissonNL::functionspace_type::mesh_support_vector_type
PoissonNL::functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
{
    auto elecDomain = markedelements(mesh, M_materials.markersWithPhysic("electric"));
    auto therDomain = markedelements(mesh, M_materials.markersWithPhysic("thermic"));
    auto suppElec = std::make_shared<MeshSupport<mesh_type>>(mesh,elecDomain);
    auto suppTher = std::make_shared<MeshSupport<mesh_type>>(mesh,therDomain);
    return fusion::make_vector(suppElec,suppTher);
}

void PoissonNL::initModel()
{
    auto propertyPath = Environment::expand( soption("poisson.filename"));
    this->addModelFile("property-file", propertyPath);
    M_modelProps = std::make_shared<prop_type>(propertyPath);

    auto parameters = M_modelProps->parameters();
    int nbCrbParameters = count_if(parameters.begin(), parameters.end(), [] (auto const& p)
                                   {
                                       return p.second.hasMinMax();
                                   });
    Dmu->setDimension(nbCrbParameters);
    auto mu_min = Dmu->element();
    auto mu_max = Dmu->element();
    int i = 0;
    for( auto const& [key,parameter] : parameters )
    {
        if( parameter.hasMinMax() )
        {
            mu_min(i) = parameter.min();
            mu_max(i) = parameter.max();
            Feel::cout << "Found parameter " << key
                       << " with min/max : " << mu_min(i) << "/" << mu_max(i) << std::endl;
            Dmu->setParameterName(i++, key );
        }
    }
    Dmu->setMin(mu_min);
    Dmu->setMax(mu_max);
    M_mu = Dmu->element();

    M_materials = M_modelProps->materials().materialWithPhysic(std::vector<std::string>({"electric","thermic"}));
    M_elecMaterials = M_modelProps->materials().materialWithPhysic("electric");
    M_therMaterials = M_modelProps->materials().materialWithPhysic("thermic");
    M_nbElecMat = M_elecMaterials.size();
    M_nbTherMat = M_therMaterials.size();
    auto bdConditions = M_modelProps->boundaryConditions2();
    auto potentialDirichlet = bdConditions.boundaryConditions("potential","Dirichlet");
    auto temperatureRobin = bdConditions.boundaryConditions("temperature","Robin");
    M_nbPotDir = potentialDirichlet.size();
    M_nbTempRobin = temperatureRobin.size();

    if( !M_mesh )
        M_mesh = loadMesh( new mesh_type );
    this->setFunctionSpaces(functionspace_type::New( _mesh=M_mesh, _range=this->functionspaceMeshSupport( M_mesh ) ) );
    Feel::cout << "Potential nDof  : " << Xh->template functionSpace<0>()->nDof() << std::endl
               << "Temperature nDof: " << Xh->template functionSpace<1>()->nDof() << std::endl;

    M_u = Xh->element();

    auto mu = this->Dmu->min();
    M_InitialGuess.resize(2);
    this->M_betaInitialGuess.resize(2);
    for( int q = 0; q < 2; ++q )
    {
        M_InitialGuess[q].resize(1, Xh->elementPtr());
        this->M_betaInitialGuess[q].resize(1);
    }
    auto VTInit = this->solveLinear(mu);
    M_InitialGuess[0][0]->template element<0>() = VTInit.template element<0>();
    M_InitialGuess[1][0]->template element<1>() = VTInit.template element<1>();

    auto Pset = this->Dmu->sampling();

    int trainsetEimSize = ioption("poisson.trainset-eim-size");
    std::string supersamplingname = (boost::format("DmuEim-Ne%1%-generated-by-master-proc")
                                     % trainsetEimSize ).str();
    std::ifstream file ( supersamplingname );
    bool all_proc_same_sampling=true;
    if( ! file )
    {
        Pset->randomize( trainsetEimSize , all_proc_same_sampling , supersamplingname );
        Pset->writeOnFile( supersamplingname );
    }
    else
    {
        Pset->clear();
        Pset->readFromFile(supersamplingname);
    }

    // eim
    auto T0 = cst(293.0);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto name = "eim_k_"+key;
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto k = sigma*L*_e1;
        auto Eh = eim_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.meshMarkers()) );
        // auto Eh = eim_space_type::New( _mesh=M_mesh );

        auto eim_k = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                          _element=M_u.template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=Eh,
                          _name=name,
                          _sampling=Pset );
        this->addEim( eim_k );
        Feel::cout << tc::green << name << " dimension: " << eim_k->mMax() << tc::reset << std::endl;
    }

    for( auto const& [key,mat] : M_elecMaterials )
    {
        auto name = "eim_sigma_"+key;
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto Eh = eim_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, mat.meshMarkers()) );
        // auto Eh = eim_space_type::New( _mesh=M_mesh );
        auto eim_sigma = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                              _element=M_u.template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma,
                              _space=Eh,
                              _name=name,
                              _sampling=Pset );
        this->addEim( eim_sigma );
        Feel::cout << tc::green << name << " dimension: " << eim_sigma->mMax() << tc::reset << std::endl;
    }

    auto gradgrad = _e2v*trans(_e2v);
    auto Jh = eimd_space_type::New( _mesh=M_mesh, _range=markedelements(M_mesh, M_materials.markersWithPhysic("electric")) );
    auto eim_grad = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                         _element=M_u.template element<0>(),
                         _parameter=M_mu,
                         _expr=gradgrad,
                         _space=Jh,
                         _name="eim_grad",
                         _sampling=Pset );
    this->addEimDiscontinuous( eim_grad );
    Feel::cout << tc::green << "eim_grad dimension: " << eim_grad->mMax() << tc::reset << std::endl;

    auto U = Xh->element();
    auto u1 = U.template element<0>();
    auto u2 = U.template element<1>();

    // Energy matrix
    auto m = form2(_test=Xh, _trial=Xh);
    m = integrate( elements(support(Xh->template functionSpace<0>())),
                   inner(gradt(u1),grad(u1)) );
    m += integrate( elements(support(Xh->template functionSpace<1>())),
                    inner(gradt(u2),grad(u2)) );
    M_energy_matrix = m.matrixPtr();

    auto mm = form2(_test=Xh, _trial=Xh);
    mm = integrate( elements(support(Xh->template functionSpace<0>())),
                    inner(idt(u1),id(u1)) );
    mm += integrate( elements(support(Xh->template functionSpace<1>())),
                     inner(idt(u2),id(u2)) );
    M_mass_matrix = mm.matrixPtr();

    this->assemble();
}

void PoissonNL::setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir )
{
    auto const& ptreeEim = ptree.get_child( "eim" );
    this->scalarContinuousEim().clear();
    this->scalarDiscontinuousEim().clear();
    auto T0 = cst(293.0);
    for( auto const& [key,mat] : M_therMaterials )
    {
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto k = sigma*L*_e1;
        eim_space_ptrtype Eh;
        auto name = "eim_k_"+key;
        auto const& ptreeEimK = ptreeEim.get_child( name );
        std::string dbNameEimK = ptreeEimK.template get<std::string>( "database-filename" );
        auto eim_k = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                          _element=M_u.template element<1>(),
                          _parameter=M_mu,
                          _expr=k,
                          _space=Eh,
                          _name=name,
                          _filename=dbNameEimK,
                          _directory=dbDir );
        this->addEim( eim_k );
        Feel::cout << tc::green << name << " dimension: " << eim_k->mMax() << tc::reset << std::endl;
    }
    for( auto const& [key,mat] : M_elecMaterials )
    {
        eim_space_ptrtype Eh;
        auto name = "eim_sigma_"+key;
        auto sigma0 = cst_ref(M_mu.parameterNamed(mat.getString("misc.sigmaKey")));
        auto alpha = cst_ref(M_mu.parameterNamed(mat.getString("misc.alphaKey")));
        auto L = cst_ref(M_mu.parameterNamed("L"));
        auto sigma = sigma0/(cst(1.) + alpha*(_e1-T0));
        auto const& ptreeEimSigma = ptreeEim.get_child( name );
        std::string dbNameEimSigma = ptreeEimSigma.template get<std::string>( "database-filename" );
        auto eim_sigma = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                              _element=M_u.template element<1>(),
                              _parameter=M_mu,
                              _expr=sigma,
                              _space=Eh,
                              _name=name,
                              _filename=dbNameEimSigma,
                              _directory=dbDir );
        this->addEim( eim_sigma );
        Feel::cout << tc::green << name << " dimension: " << eim_sigma->mMax() << tc::reset << std::endl;
    }

    auto const& ptreeEimGrad = ptreeEim.get_child( "eim_grad" );
    std::string dbNameEimGrad = ptreeEimGrad.template get<std::string>( "database-filename" );
    auto gradgrad = _e2v*trans(_e2v);
    eimd_space_ptrtype Jh;
    auto eim_grad = eim( _model=std::dynamic_pointer_cast<PoissonNL>(this->shared_from_this() ),
                         _element=M_u.template element<0>(),
                         _parameter=M_mu,
                         _expr=gradgrad,
                         _space=Jh,
                         _name="eim_grad",
                         _filename=dbNameEimGrad,
                         _directory=dbDir );
    this->addEimDiscontinuous( eim_grad );
    Feel::cout << tc::green << "eim_grad dimension: " << eim_grad->mMax() << tc::reset << std::endl;

    this->resize();
}

double PoissonNL::output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve)
{
    auto mesh = Xh->mesh();
    double output=0;
    if ( output_index == 0 )
    {
        for ( int q = 0; q < Ql(0); q++ )
        {
            element_ptrtype eltF( new element_type( Xh ) );
            *eltF = *M_Fq[output_index][q];
            output += M_betaFq[output_index][q]*dot( *eltF, u );
            //output += M_betaFqm[output_index][q][m]*dot( M_Fqm[output_index][q][m], U );
        }
    }
    else
        throw std::logic_error( "[Heat2d::output] error with output_index : only 0 or 1 " );
    return output;
}

// FEELPP_CRB_PLUGIN( PoissonNL, poissonnl )
}
