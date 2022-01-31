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

#ifndef BIOTSAVART_RB_HPP
#define BIOTSAVART_RB_HPP 1

#include <feel/feelmor/crbdb.hpp>
#include <feel/feelmor/crbdata.hpp>
#include <feel/feelmor/crbelementsdb.hpp>
#include <feel/feelcore/core.hpp>
#include <feel/feeldiscr/expansion.hpp>

namespace Feel
{
/**
 * \class CRBPluginBase
 */
template<typename ModelT>
class CRBPluginBase : public CRBDB
{
public:
    using self_type = CRBPluginBase;
    using self_ptrtype = std::shared_ptr<self_type>;

    using parameterspace_type = ParameterSpaceX;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameter_type = parameterspace_type::element_type;

    using model_type = ModelT;
    using model_ptrtype = std::shared_ptr<model_type>;

    CRBPluginBase() = delete;
    CRBPluginBase( std::string const& name, std::string const& ext,
                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : CRBDB( name, ext, worldComm )
        {
            Feel::cout << tc::red << "CRBPluginBase constructor with uuid: " << id() << tc::reset << std::endl;
        }
    CRBPluginBase( std::string const& name, std::string const& ext, uuids::uuid const& i,
                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : CRBDB( name, ext, i, worldComm )
        {
            Feel::cout << tc::red << "CRBPluginBase constructor with uuid: " << id() << tc::reset << std::endl;
        }
    CRBPluginBase( CRBPluginBase const& ) = default;
    CRBPluginBase& operator=( CRBPluginBase const& o ) = default;

    virtual parameterspace_ptrtype Dmu() const = 0;
    virtual CRBResults run( parameter_type const& mu, vectorN_type& time, double eps, int N, bool print_rb_matrix ) = 0;
    virtual std::vector<CRBResults> run( std::vector<parameter_type> const& mu, double eps, int N, bool print_rb_matrix ) = 0;
    virtual model_ptrtype model() const = 0;
};

template <typename BiotSavartModelT>
class BiotSavartRB : public CRBPluginBase<BiotSavartModelT>
{
public:
    using model_type = BiotSavartModelT;
    using model_ptrtype = std::shared_ptr<model_type>;
    using super = CRBPluginBase<model_type>;
    using self_type = BiotSavartRB;
    using self_ptrtype = std::shared_ptr<self_type>;

    using parameterspace_type = ParameterSpaceX;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameter_type = parameterspace_type::element_type;

    using mesh_ptrtype = typename model_type::mesh_ptrtype;
    using element_type = typename model_type::element_type;
    using vt_element_type = typename model_type::vt_element_type;

    using crb_elements_db_type = CRBElementsDB<model_type>;
    using crb_elements_db_ptrtype = std::shared_ptr<crb_elements_db_type>;

private:
    model_ptrtype M_model;
    crb::stage M_stage;

    int M_N;
    crb_elements_db_type M_elements_db;
    bool M_rebuild;

protected:
    BiotSavartRB( std::string const& name,
                  model_ptrtype const& model,
                  crb::stage = crb::stage::online,
                  std::string const& prefixExt = "",
                  std::string const& prefix = "biotsavart");

public:
    static self_ptrtype New( std::string const& name, crb::stage stage = crb::stage::online, std::string const& prefix="biotsavart" );
    static self_ptrtype New( std::string const& name,
                             model_ptrtype const& model,
                             crb::stage stage = crb::stage::online,
                             std::string const& prefixExt = "",
                             std::string const& prefix = "biotsavart");
    void init();

    void saveRB();
    void saveDB() override;
    void saveJson();
    void setupOfflineFromDB();
    void loadDB( std::string const& filename, crb::load l ) override;
    void loadJson( std::string const& filename );

    parameterspace_ptrtype Dmu() const override { return M_model->parameterSpace(); }
    CRBResults run( parameter_type const& mu, vectorN_type& time, double eps=1e-6, int N=-1, bool print_rb_matrix=false ) override;
    std::vector<CRBResults> run( std::vector<parameter_type> const& mu, double eps=1e-6, int N=-1, bool print_rb_matrix=false ) override;
    model_ptrtype model() const override { return M_model; }
    element_type expansion( vectorN_type const& b, int N = -1 ) const;
    vt_element_type expansionVT( vectorN_type const& vtN, int N = -1 ) const;

    void exportBasis() const;

    void offline();
    vectorN_type onlineVT( parameter_type const& mu );
    vectorN_type online( parameter_type const& mu );
    vectorN_type online( parameter_type const& mu, vectorN_type const& vtN );

    mesh_ptrtype const& mesh() const { return M_model->mesh(); }
    parameter_type paramFromProperties() const { return M_model->paramFromProperties(); }
};

template <typename BiotSavartModelT>
typename BiotSavartRB<BiotSavartModelT>::self_ptrtype
BiotSavartRB<BiotSavartModelT>::New( std::string const& name, crb::stage stage, std::string const& prefix )
{
    return New( name, model_type::New(stage, prefix), stage );
}

template <typename BiotSavartModelT>
typename BiotSavartRB<BiotSavartModelT>::self_ptrtype
BiotSavartRB<BiotSavartModelT>::New( std::string const& name,
                                     model_ptrtype const& model,
                                     crb::stage stage,
                                     std::string const& prefixExt,
                                     std::string const& prefix)
{
    auto bs = std::shared_ptr<self_type>( new self_type(name, model, stage, prefixExt, prefix) );
    if( stage == crb::stage::offline )
        bs->init();
    return bs;
}

template <typename BiotSavartModelT>
BiotSavartRB<BiotSavartModelT>::BiotSavartRB( std::string const& name,
                                              model_ptrtype const& model,
                                              crb::stage stage,
                                              std::string const& prefixExt,
                                              std::string const& prefix)
    : super( name, prefixvm(prefixExt, "biotsavart"), model->uuid()),
      M_model(model),
      M_stage(stage),
      M_N(0),
      M_elements_db( name, prefixvm(prefixExt, "bs-elements")),
      M_rebuild( boption("biotsavart.rebuild-database") )
{
    Feel::cout << tc::red << "BiotSavartRB constructor with uuid: " << this->id() << tc::reset << std::endl;
    M_elements_db.setModel( model );
    if( M_stage == crb::stage::online )
        this->loadDB( (this->dbLocalPath()/fs::path(this->jsonFilename())).string(), crb::load::all );
}

template <typename BiotSavartModelT>
void BiotSavartRB<BiotSavartModelT>::init()
{
    if( !M_rebuild && !fs::exists(this->dbLocalPath()/fs::path(this->jsonFilename())) )
    {
        M_rebuild = true;
        this->worldComm().barrier();
    }
    if( !M_rebuild )
    {
        this->setupOfflineFromDB();
    }
    else
    {
        M_elements_db.setId( this->id() );
    }
    Feel::cout << "Use DB id " << this->id() << std::endl;
    if ( M_N == 0 )
    {
        Feel::cout << "Databases does not exist or incomplete -> Start from the begining" << std::endl;;
        LOG( INFO ) << "Databases does not exist or incomplete -> Start from the begining";
    }
}

template <typename BiotSavartModelT>
void
BiotSavartRB<BiotSavartModelT>::saveRB()
{
    M_elements_db.setWn( std::make_tuple( M_model->rBFunctionSpace()->primalRB() , M_model->rBFunctionSpace()->dualRB() ) );

    M_elements_db.saveDB();
}

template <typename BiotSavartModelT>
void
BiotSavartRB<BiotSavartModelT>::saveDB()
{
    saveJson();
}

template <typename BiotSavartModelT>
void
BiotSavartRB<BiotSavartModelT>::saveJson()
{
    if ( this->worldComm().isMasterRank() )
    {
        std::string filenameJson = (this->dbLocalPath()/fs::path(this->jsonFilename())).string();
        std::cout << "saveDB: " << filenameJson << std::endl;

        boost::property_tree::ptree ptree;
        ptree.add( "uuid", this->idStr() );

        boost::property_tree::ptree ptreeReducedBasisSpace;
        std::string meshFilename = (boost::format("%1%_mesh_p%2%.json")%this->name() %this->worldComm().size()).str();
        ptreeReducedBasisSpace.add( "mesh-filename",meshFilename );
        ptreeReducedBasisSpace.add( "database-filename", M_elements_db.dbFilename() );
        ptreeReducedBasisSpace.add( "dimension", M_N );
        if ( M_model && M_model->rBFunctionSpace() && M_model->rBFunctionSpace()->functionSpace() )
        {
            ptreeReducedBasisSpace.add( "mesh-context",M_model->rBFunctionSpace()->functionSpace()->mesh()->components().context() );

            auto feSpace = M_model->rBFunctionSpace()->functionSpace();
            boost::property_tree::ptree ptreeFiniteElementSpace;
            ptreeFiniteElementSpace.add( "dimension", feSpace->nDof() );
            ptreeFiniteElementSpace.add( "basis-name", feSpace->basisName() );
            ptreeReducedBasisSpace.add_child( "finite-element-space", ptreeFiniteElementSpace );
        }
        ptree.add_child( "reduced-basis-space", ptreeReducedBasisSpace );

        boost::property_tree::ptree ptreeCrb;//Database;
        ptreeCrb.add( "dimension", M_N );
        ptreeCrb.add( "name", this->name() );
        ptree.add_child( "crb", ptreeCrb );

        write_json( filenameJson, ptree );
    }
}

template <typename BiotSavartModelT>
void BiotSavartRB<BiotSavartModelT>::setupOfflineFromDB()
{
    this->loadDB( (this->dbLocalPath()/fs::path(this->jsonFilename())).string(), crb::load::all );
    if( M_elements_db.loadDB() )
    {
        Feel::cout << "Database for basis functions " << M_elements_db.lookForDB()
                   << " available and loaded with " << M_N << " basis" << std::endl;
        LOG(INFO) << "Database for basis functions " << M_elements_db.lookForDB()
                  << " available and loaded with " << M_N << " basis";
        auto basis_functions = M_elements_db.wn();
        M_model->rBFunctionSpace()->setBasis( basis_functions );
    }
}

template <typename BiotSavartModelT>
void BiotSavartRB<BiotSavartModelT>::loadDB( std::string const& filename, crb::load l )
{
    auto fname = this->db( filename );
    this->loadJson( fname.string() );
}

template <typename BiotSavartModelT>
void BiotSavartRB<BiotSavartModelT>::loadJson( std::string const& filename )
{
    std::string dbDir = fs::path( filename ).parent_path().string();

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    boost::property_tree::ptree ptree;
    std::istringstream istr( json_str_wo_comments );
    boost::property_tree::read_json( istr, ptree );

    auto i = ptree.template get<std::string>( "uuid" );
    this->setId( boost::lexical_cast<uuids::uuid>( i ) );

    auto const& ptreeCrb = ptree.get_child( "crb" );
    M_N = ptreeCrb.template get<int>( "dimension" );
    M_model->setN( M_N );
    this->setName( ptreeCrb.template get<std::string>( "name" ) );

    auto const& ptreeReducedBasisSpace = ptree.get_child( "reduced-basis-space" );
    M_elements_db.setup( ptreeReducedBasisSpace, dbDir );
}

template <typename BiotSavartModelT>
void
BiotSavartRB<BiotSavartModelT>::offline()
{
    M_model->crbOffline();
    // try to reload elements
    if( M_rebuild )
    {
        Feel::cout << "Database for biotsavart not found, computing the database" << std::endl;
        for( int n = 0; n < M_model->Qb(); ++n )
        {
            for( int m = 0; m < M_model->mMaxB(); ++m )
            {
                M_N = M_model->computeIntegrals( n, m );
                this->saveDB();
                this->saveRB();
            }
        }
    }
    Feel::cout << tc::green << "BiotSavart dimension: " << M_N << tc::reset << std::endl;
}

template <typename BiotSavartModelT>
vectorN_type
BiotSavartRB<BiotSavartModelT>::onlineVT( parameter_type const& mu )
{
    return M_model->computeCoeffVT( mu );
}

template <typename BiotSavartModelT>
vectorN_type
BiotSavartRB<BiotSavartModelT>::online( parameter_type const& mu )
{
    return M_model->computeCoeff( mu );
}

template <typename BiotSavartModelT>
vectorN_type
BiotSavartRB<BiotSavartModelT>::online( parameter_type const& mu, vectorN_type const& vtN )
{
    return M_model->computeCoeff( mu, vtN );
}

template <typename BiotSavartModelT>
CRBResults
BiotSavartRB<BiotSavartModelT>::run( parameter_type const& mu, vectorN_type& time, double eps, int N, bool print_rb_matrix )
{
    auto res = CRBResults();
    return res;
}

template <typename BiotSavartModelT>
std::vector<CRBResults>
BiotSavartRB<BiotSavartModelT>::run( std::vector<parameter_type> const& mu, double eps, int N, bool print_rb_matrix )
{
    auto res = std::vector<CRBResults>(mu.size());
    return res;
}

template <typename BiotSavartModelT>
typename BiotSavartRB<BiotSavartModelT>::element_type
BiotSavartRB<BiotSavartModelT>::expansion( vectorN_type const& b, int N ) const
{
    auto WN = M_model->rBFunctionSpace()->primalRB();
    int Nwn;
    if( N == -1 )
        Nwn = std::min<int>( WN.size(), b.size() );
    else
        Nwn = std::min<int>( { static_cast<int>(WN.size()), static_cast<int>(b.size()), N } );
    return Feel::expansion( WN, b, Nwn );
}

template <typename BiotSavartModelT>
typename BiotSavartRB<BiotSavartModelT>::vt_element_type
BiotSavartRB<BiotSavartModelT>::expansionVT( vectorN_type const& vtN, int N ) const
{
    return M_model->expansionVT( vtN, N );
}

template <typename BiotSavartModelT>
void
BiotSavartRB<BiotSavartModelT>::exportBasis() const
{
    auto WN = M_model->rBFunctionSpace()->primalRB();
    int N = WN.size();
    auto e = exporter( _mesh=M_model->mesh(), _name="basis" );
    for( int n = 0; n < N; ++n )
        e->add( (boost::format("basis_%1%") % n ).str(), *WN[n] );
    e->save();
}



}

#endif
