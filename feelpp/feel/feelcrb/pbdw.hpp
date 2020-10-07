/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Romain Hild <romain.hild@cemosis.fr>
       Date: 2020-09-30

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef _FEELPP_GEIM_HPP
#define _FEELPP_GEIM_HPP 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcore/unwrapptr.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feelcore/json.hpp>

#include <limits>

namespace Feel
{

/**
 * Class to support Parameterized Background Data-Weak method
 */
template<typename RBSpace>
class PBDW : CRBDB
{
public:
    using super_type = CRBDB;
    using reducedspace_type = RBSpace;
    using reducedspace_ptrtype = std::shared_ptr<reducedspace_type>;
    using space_type = typename reducedspace_type::fespace_type;
    using space_ptrtype = std::shared_ptr<space_type>;
    using mesh_type = typename space_type::mesh_type;
    using node_t = typename mesh_type::node_type;
    using element_type = typename space_type::element_type;
    using vectorN_type = Eigen::VectorXd;
    using matrixN_type = Eigen::MatrixXd;
    using sensorbase_type = SensorBase<space_type>;
    using sensorbase_ptrtype = std::shared_ptr<sensorbase_type>;
    static const int nDim = space_type::mesh_type::nDim;

    /**
     * Constructor for the online phase
     * @param name Name of pbdw
     * @param l Loading type (rb,fe,all)
     * @param uuid Uuid to use for the db
     */
    explicit PBDW(std::string const& name,
                  crb::load l = crb::load::rb,
                  uuids::uuid const& uuid = uuids::nil_uuid());
    /**
     * Constructor for the offline phase
     * @param name Name of pbdw
     * @param XR Reduced Basis
     * @param sigmas Sensors to use
     * @param uuid Uuid to use for the db
     */
    PBDW(std::string const& name,
         reducedspace_ptrtype const& XR,
         std::vector<sensorbase_ptrtype> const& sigmas,
         uuids::uuid const& uuid = uuids::nil_uuid());
    int dimensionN() const { return M_N; } /**< Dimension of Reduced Basis */
    int dimensionM() const { return M_M; } /**< Number of sensors */
    int dimension() const { return M_N+M_M; } /**< Dimension of PBDW */
    space_ptrtype functionSpace() const { return M_XR->functionSpace(); } /**< Function Space */
    matrixN_type matrix() const { return M_matrix; } /**< Matrix of PBDW */
    void offline(); /** Do offline phase */
    /**
     * Do online phase
     * @param yobs Observation of sensors
     * @return Coefficients of solution
     */
    vectorN_type online(vectorN_type const& yobs) const; /**< Do online phase */
    /**
     * Retrieve solution
     * @param yobs Observation of sensors
     * @return Solution
     */
    element_type solution(vectorN_type const& yobs) const; /**< Retrieve solution */

protected:
    void loadDB( std::string const& filename, crb::load l ) override;
    void saveDB() override;


private:
    std::string M_name;
    reducedspace_ptrtype M_XR;
    std::vector<sensorbase_ptrtype> M_sigmas;
    matrixN_type M_matrix;
    int M_M;
    int M_N;

    bool M_rebuildDb;
    int M_dbLoad;
    std::string M_dbFilename;
    std::string M_dbId;
    crb::stage M_stage;
};

template<typename RBSpace>
PBDW<RBSpace>::PBDW(std::string const& name,
                    crb::load l,
                    uuids::uuid const& uuid):
    super_type(name, "pbdw", uuid),
    M_name(name),
    M_rebuildDb(boption("pbdw.rebuild-database")),
    M_dbLoad(ioption("pbdw.db.load")),
    M_dbFilename(soption("pbdw.db.filename")),
    M_dbId(soption("pbdw.db.id")),
    M_stage(crb::stage::online)
{
    if( ! this->findDBUuid(M_dbLoad, M_dbLoad ? M_dbId : M_dbFilename) )
        throw std::invalid_argument("Database not found during online phase");

    this->loadDB(this->absoluteDbFilename(), l );
}

template<typename RBSpace>
PBDW<RBSpace>::PBDW(std::string const& name,
                    reducedspace_ptrtype const& XR,
                    std::vector<sensorbase_ptrtype> sigmas,
                    uuids::uuid const& uuid):
    super_type(name, "pbdw", uuid),
    M_name(name),
    M_XR(XR),
    M_sigmas(sigmas),
    M_M(M_sigmas.size()),
    M_N(M_XR->size()),
    M_rebuildDb(boption("pbdw.rebuild-database")),
    M_dbLoad(ioption("pbdw.db.load")),
    M_dbFilename(soption("pbdw.db.filename")),
    M_dbId(soption("pbdw.db.id")),
    M_stage(crb::stage::offline)
{
    if( ! this->findDBUuid(M_dbLoad, M_dbLoad ? M_dbId : M_dbFilename) )
        this->setDBDirectory(Environment::randomUUID(true));

    if( ! M_rebuildDb )
        this->loadDB(this->absoluteDbFilename(), crb::load::all );
}

template<typename RBSpace>
void
PBDW<RBSpace>::offline()
{
    Feel::cout << "offline phase start with M=" << M_M << " and N=" << M_N << std::endl;
    M_matrix = matrixN_type::Zero(M_M+M_N, M_M+M_N);
    for(int i = 0; i < M_M; ++i )
    {
        for(int j = 0; j < i; ++j )
        {
            M_matrix(i, j) = inner_product(M_sigmas[i]->containerPtr(), M_sigmas[j]->containerPtr());
            M_matrix(j, i) = M_matrix(i, j);
        }
        M_matrix(i, i) = inner_product(M_sigmas[i]->containerPtr(), M_sigmas[i]->containerPtr());
        for(int j = 0; j < M_N; ++j )
            M_matrix(i, M_M+j) = (*M_sigmas[i])(M_XR->primalBasisElement(j));
    }
    M_matrix.bottomLeftCorner(M_N, M_M) = M_matrix.topRightCorner(M_M, M_N).transpose();

    this->saveDB();
}

template<typename RBSpace>
typename PBDW<RBSpace>::vectorN_type
PBDW<RBSpace>::online(vectorN_type const& yobs) const
{
    vectorN_type yobs2 = vectorN_type::Zero(this->dimension());
    yobs2.head(this->dimensionM()) = yobs;
    vectorN_type vn = M_matrix.colPivHouseholderQr().solve(yobs2);
    return vn;
}

template<typename RBSpace>
typename PBDW<RBSpace>::element_type
PBDW<RBSpace>::solution(vectorN_type const& yobs) const
{
    int M = this->dimensionM();
    int N = this->dimensionN();
    auto vn = this->online(yobs);
    auto wn = M_XR->primalRB();
    auto Xh = M_XR->functionSpace();
    auto I = Xh->element();
    I.zero();
    for(int i = 0; i < M; ++i )
        I.add(vn(i), Xh->element(M_sigmas[i]->containerPtr()));
    for(int i = 0; i < N; ++i )
        I.add(vn(M+i), unwrap_ptr(wn[i]));
    return I;
}

template<typename RBSpace>
void
PBDW<RBSpace>::saveDB()
{
    if ( this->worldComm().isMasterRank() )
    {
        if( fs::exists(this->dbLocalPath()) )
        {
            if( !fs::is_directory(this->dbLocalPath()) )
                throw std::exception();
        }
        else
            fs::create_directories(this->dbLocalPath());

        json j;
        j["name"] = this->name();
        std::time_t t = std::time(nullptr);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&t), "%c %Z");
        j["date"] = ss.str();
        j["mesh"] = this->absoluteMeshFilename();
        std::ofstream o(this->absoluteJsonFilename());
        o << j << std::endl;

        fs::ofstream ofs( this->absoluteDbFilename() );
        if ( ofs )
        {
            Feel::cout << "saving DB at " << this->absoluteDbFilename() << std::endl;
            boost::archive::binary_oarchive oa( ofs );
            oa << this->dimensionN();
            oa << this->dimensionM();
            oa << M_matrix;
        }
    }
    if( ! fs::exists(this->absoluteMeshFilename()) )
    {
        Feel::cout << "Saving mesh to " << this->absoluteMeshFilename() << std::endl;
        M_XR->functionSpace()->mesh()->saveHDF5( this->absoluteMeshFilename() );
    }
    fs::ofstream ofsp( this->absoluteDbFilenameProc() );
    if( ofsp )
    {
        boost::archive::binary_oarchive oa( ofsp );
        for( int n = 0; n < M_N; ++n )
            oa << M_XR->primalBasisElement(n);
        for( int m = 0; m < M_M; ++m )
        {
            auto t = M_sigmas[m]->type();
            oa << t;
            oa << M_sigmas[m]->name();
            oa << M_sigmas[m]->position();
            if( t == "pointwise" )
            {
            }
            else if( t == "gaussian" )
            {
                auto s = std::dynamic_pointer_cast<SensorGaussian<space_type>>(M_sigmas[m]);
                oa << s->radius();
            }
        }
    }
}

template<typename RBSpace>
void
PBDW<RBSpace>::loadDB( std::string const& filename, crb::load l )
{
    if( !fs::exists(filename) )
        return;

    fs::ifstream ifs( filename );
    if ( ifs )
    {
        Feel::cout << "loading DB at " << filename << std::endl;
        boost::archive::binary_iarchive ia( ifs );
        if( l > crb::load::none )
        {
            ia >> M_N;
            ia >> M_M;
            M_matrix.resize(M_M+M_N, M_M+M_N);
            ia >> M_matrix;
        }
    }
    if( l > crb::load::rb )
    {
        if( ! M_XR )
        {
            std::ifstream i(this->absoluteJsonFilename());
            json j = json::parse(i);
            auto meshfilename = j["mesh"].get<std::string>();
            auto mesh = loadMesh(_mesh=new mesh_type, _filename=meshfilename);
            auto Xh = space_type::New(mesh);
            M_XR = std::make_shared<reducedspace_type>(Xh);
            fs::ifstream ifsp( this->absoluteDbFilenameProc() );
            if( ifsp )
            {
                boost::archive::binary_iarchive ia( ifsp );
                auto u = Xh->element();
                for( int n = 0; n < M_N; ++n )
                {
                    ia >> u;
                    M_XR->addPrimalBasisElement(u);
                }
                node_t n(nDim);
                for( int m = 0; m < M_M; ++m )
                {
                    std::string type;
                    ia >> type;
                    std::string name;
                    ia >> name;
                    node_t n(nDim);
                    ia >> n;
                    if( type == "pointwise" )
                    {
                        auto s = std::make_shared<SensorPointwise<space_type> >(Xh, n, name);
                        M_sigmas.push_back(s);
                    }
                    else if( type == "gaussian" )
                    {
                        double r;
                        ia >> r;
                        auto s = std::make_shared<SensorGaussian<space_type> >(Xh, n, r, name);
                        M_sigmas.push_back(s);
                    }
                }
            }
        }
    }
}

} // namespace Feel

#endif
