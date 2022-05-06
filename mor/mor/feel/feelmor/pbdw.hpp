/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Romain Hild <romain.hild@cemosis.fr>
       Date: 2020-09-30

  Copyright (C) 2020 Feel++ Consortium
  Copyright (C) 2020 University of Strasbourg

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
#ifndef _FEELPP_PBDW_HPP
#define _FEELPP_PBDW_HPP 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelmor/crbdb.hpp>
#include <feel/feelcore/unwrapptr.hpp>
#include <feel/feelmor/reducedbasisspace.hpp>
#include <feel/feeldiscr/sensors.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelfilters/loadmesh.hpp>

#include <limits>

namespace Feel
{

/**
 * Class that implements the online phase of the Parameterized Background Data-Weak method.
 */
class PBDWOnline : public CRBDB
{
public:
    using super_type = CRBDB;
    using vectorN_type = Eigen::VectorXd;
    using matrixN_type = Eigen::MatrixXd;

    /**
     * Constructor for the online phase from PBDW
     * @param name Name of pbdw
     * @param uuid Uuid to use for the db
     */
    explicit PBDWOnline(std::string const& name, uuids::uuid const& uuid);
    /**
     * Constructor for the online phase
     * @param name Name of pbdw
     * @param dbLoad Loading type for the DB
     * @param dbFilename Filename of the DB for load type filename
     * @param dbIf If of the DB for load type id
     */
    explicit PBDWOnline(std::string const& name,
                        int dbLoad = ioption("pbdw.db.load"),
                        std::string const& dbFilename = soption("pbdw.db.filename"),
                        std::string const& dbId = soption("pbdw.db.id"));
    int dimensionN() const { return M_N; } /**< Dimension of Reduced Basis */
    int dimensionM() const { return M_M; } /**< Number of sensors */
    int dimension() const { return M_N+M_M; } /**< Dimension of PBDW */
    int dimensionF() const { return M_Nl; } /**< Dimension of outputs */
    matrixN_type matrix() const { return M_matrix; } /**< Matrix of PBDW */
    matrixN_type matrixF() const { return M_F; } /**< Matrix of outputs */
    std::vector<std::string> sensorNames() const { return M_sensorNames; } /**< Names of the sensors */
    std::vector<std::string> outputNames() const { return M_outputNames; } /**< Names of the outputs */
    /**
     * Do online phase
     * @param yobs Observation of sensors
     * @return Coefficients of solution
     */
    matrixN_type online(matrixN_type const& yobs) const;
    /**
     * Get outputs
     * @param yobs Observation of sensors
     * @return Outputs
     */
    matrixN_type outputs(matrixN_type const& yobs) const;
    /**
     * Do online phase
     * @param yobs Observation of sensors
     * @param sensors Index of sensors to use
     * @param toComplete complete indexes to M+N dimension
     * @return Coefficients of solution
     */
    matrixN_type online(matrixN_type const& yobs, std::vector<int> const& sensors, bool toComplete = true) const;
    /**
     * Get outputs
     * @param yobs Observation of sensors
     * @param sensors Index of sensors to use
     * @param toComplete complete indexes to M+N dimension
     * @return Outputs
     */
    matrixN_type outputs(matrixN_type const& yobs, std::vector<int> const& sensors, bool toComplete = true) const;
    /**
     * Do online phase
     * @param yobs Observation of sensors
     * @param sensors Names of sensors to use
     * @return Coefficients of solution
     */
    matrixN_type online(matrixN_type const& yobs, std::vector<std::string> const& sensors) const;
    /**
     * Get outputs
     * @param yobs Observation of sensors
     * @param sensors Names of sensors to use
     * @return Outputs
     */
    matrixN_type outputs(matrixN_type const& yobs, std::vector<std::string> const& sensors) const;

protected:
    void loadDB( std::string const& filename, crb::load l ) override;
    std::vector<int> idsWith( std::vector<int> const& sensors ) const;
    std::vector<int> idsWith( std::vector<std::string> const& sensors ) const;

    std::string M_name;
    matrixN_type M_matrix;
    matrixN_type M_F;
    int M_M;
    int M_N;
    int M_Nl;
    std::vector<std::string> M_sensorNames;
    std::vector<std::string> M_outputNames;
};

PBDWOnline::PBDWOnline(std::string const& name,
                       uuids::uuid const& uuid):
    super_type(name, "pbdw", uuid),
    M_name(name),
    M_M(0),
    M_N(0),
    M_Nl(0),
    M_matrix(matrixN_type::Zero(0,0)),
    M_F(matrixN_type::Zero(0,0))
{}

PBDWOnline::PBDWOnline(std::string const& name,
                       int dbLoad,
                       std::string const& dbFilename,
                       std::string const& dbId):
    super_type(name, "pbdw", uuids::nil_uuid()),
    M_name(name),
    M_M(0),
    M_N(0),
    M_Nl(0),
    M_matrix(matrixN_type::Zero(0,0)),
    M_F(matrixN_type::Zero(0,0))
{
    if( ! this->findDBUuid(dbLoad, dbLoad ? dbId : dbFilename) )
        throw std::invalid_argument("Database not found during online phase");

    this->loadDB(this->absoluteDbFilename(), crb::load::rb );
}

typename PBDWOnline::matrixN_type
PBDWOnline::online(matrixN_type const& yobs) const
{
    matrixN_type yobs2 = matrixN_type::Zero(this->dimension(),yobs.cols());
    yobs2.topRows(this->dimensionM()) = yobs;
    matrixN_type vn = M_matrix.colPivHouseholderQr().solve(yobs2);
    return vn;
}

typename PBDWOnline::matrixN_type
PBDWOnline::outputs(matrixN_type const& yobs) const
{
    matrixN_type coeffs = this->online(yobs);
    matrixN_type vn = M_F*coeffs;
    return vn;
}

std::vector<int>
PBDWOnline::idsWith( std::vector<int> const& sensors ) const
{
    std::vector<int> ids(sensors.size()+this->M_N);
    auto it = std::copy(sensors.begin(), sensors.end(), ids.begin());
    std::iota(it, ids.end(), this->M_M);
    return ids;
}
std::vector<int>
PBDWOnline::idsWith( std::vector<std::string> const& sensors ) const
{
    std::vector<int> ids;
    for( auto const& name : sensors )
    {
        auto it = std::find(this->M_sensorNames.begin(), this->M_sensorNames.end(), name);
        if( it != this->M_sensorNames.end() )
            ids.push_back(std::distance(this->M_sensorNames.begin(), it));
        else
            LOG(WARNING) << "sensor " << name << " does not exist !";
    }
    return idsWith(ids);
}

typename PBDWOnline::matrixN_type
PBDWOnline::online(matrixN_type const& yobs, std::vector<int> const& sensors, bool toComplete) const
{
    std::vector<int> ids;
    if( toComplete )
        ids = idsWith(sensors);
    else
        ids = sensors;
    matrixN_type yobs2 = matrixN_type::Zero(ids.size(),yobs.cols());
    yobs2.topRows(ids.size() - this->M_N) = yobs;
    matrixN_type vn = M_matrix(ids,ids).colPivHouseholderQr().solve(yobs2);
    return vn;
}

typename PBDWOnline::matrixN_type
PBDWOnline::outputs(matrixN_type const& yobs, std::vector<int> const& sensors, bool toComplete) const
{
   std::vector<int> ids;
    if( toComplete )
        ids = idsWith(sensors);
    else
        ids = sensors;
    matrixN_type coeffs = this->online(yobs, ids, false);
    matrixN_type vn = M_F(Eigen::all, ids)*coeffs;
    return vn;
}

typename PBDWOnline::matrixN_type
PBDWOnline::online(matrixN_type const& yobs, std::vector<std::string> const& sensors) const
{
    std::vector<int> ids = idsWith(sensors);
    return this->online(yobs, ids, false);
}

typename PBDWOnline::matrixN_type
PBDWOnline::outputs(matrixN_type const& yobs, std::vector<std::string> const& sensors) const
{
    std::vector<int> ids = idsWith(sensors);
    return this->outputs(yobs, ids, false);
}

void
PBDWOnline::loadDB( std::string const& filename, crb::load l )
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
            ia >> this->M_N;
            ia >> this->M_M;
            this->M_matrix.resize(this->M_M+this->M_N, this->M_M+this->M_N);
            ia >> this->M_matrix;
            ia >> this->M_Nl;
            this->M_F.resize(this->M_Nl, this->M_M+this->M_N);
            ia >> this->M_F;
            this->M_sensorNames.reserve(M_M);
            ia >> this->M_sensorNames;
            ia >> this->M_outputNames;
        }
    }
}

/**
 * Class that implements Parameterized Background Data-Weak method.
 * @tparam Reduced basis space
 */
template<typename RBSpace>
class PBDW : public PBDWOnline
{
public:
    using super_type = PBDWOnline;
    using reducedspace_type = RBSpace;
    using reducedspace_ptrtype = std::shared_ptr<reducedspace_type>;
    using space_type = typename reducedspace_type::fespace_type;
    using space_ptrtype = std::shared_ptr<space_type>;
    using element_type = typename space_type::element_type;
    using mesh_type = typename space_type::mesh_type;
    using node_t = typename mesh_type::node_type;
    using vectorN_type = Eigen::VectorXd;
    using matrixN_type = Eigen::MatrixXd;
    using sensormap_type = SensorMap<space_type>;
    using sparse_matrix_type = typename Backend<double>::sparse_matrix_type;
    using sparse_matrix_ptrtype = typename Backend<double>::sparse_matrix_ptrtype;
    using vector_ptrtype = typename Backend<double>::vector_ptrtype;
    static const int nDim = space_type::mesh_type::nDim;

    /**
     * Constructor for the online phase
     * @param name Name of pbdw
     * @param l Loading type (rb,fe,all)
     * @param uuid Uuid to use for the db
     * @param dbLoad Loading type for the DB (0: use db.filename, 1: use last DB created, 2: use last DB modified, 3: use db.id, 4: create new db)
     * @param dbFilename Filename of the DB for load type filename
     * @param dbIf If of the DB for load type id
     */
    explicit PBDW(std::string const& name,
                  crb::load l = crb::load::rb,
                  uuids::uuid const& uuid = uuids::nil_uuid(),
                  int dbLoad = ioption("pbdw.db.load"),
                  std::string const& dbFilename = soption("pbdw.db.filename"),
                  std::string const& dbId = soption("pbdw.db.id"));
    /**
     * Constructor for the offline phase
     * @param name Name of pbdw
     * @param XR Reduced Basis
     * @param sigmas Sensors to use
     * @param uuid Uuid to use for the db
     * @param rebuildDB Boolean to rebuild the database
     * @param dbLoad Loading type for the DB
     * @param dbFilename Filename of the DB for load type filename
     * @param dbIf If of the DB for load type id
     */
    PBDW(std::string const& name,
         reducedspace_ptrtype const& XR,
         sensormap_type const& sigmas,
         uuids::uuid const& uuid = uuids::nil_uuid(),
         bool rebuildDB = boption("pbdw.rebuild-database"),
         int dbLoad = ioption("pbdw.db.load"),
         std::string const& dbFilename = soption("pbdw.db.filename"),
         std::string const& dbId = soption("pbdw.db.id"));
    void setFunctionSpaceSupport(std::set<std::string> const& markers) { M_markers = markers; }
    space_ptrtype functionSpace() const { return M_XR->functionSpace(); } /**< Function Space */
    sensormap_type const& sensors() const { return M_sigmas; } /**< Sensors */
    std::vector<element_type> const& riesz() const { return M_qs; } /**< Riesz representant */
    reducedspace_ptrtype const& reducedBasis() const { return M_XR; } /**< Reduced basis */
    std::vector<vector_ptrtype> const& outputsLinear() const { return M_Fs; } /**< Linear functional for outputs */
    sparse_matrix_ptrtype initRiesz(); /**< Initialize the Riesz representation of the sensors */
    void offline(); /** Do offline phase */
    /**
     * Adds outputs to PBDW
     * @param Fs vector of functionals to apply to the basis
     */
    void setOutputs( std::vector<vector_ptrtype> const& Fs, std::vector<std::string> const& names = {} );
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
    reducedspace_ptrtype M_XR;
    sensormap_type M_sigmas;
    std::vector<element_type> M_qs;
    std::vector<vector_ptrtype> M_Fs;

    bool M_rebuildDb;
    crb::stage M_stage;
    std::set<std::string> M_markers;
};

template<typename RBSpace>
PBDW<RBSpace>::PBDW(std::string const& name,
                    crb::load l,
                    uuids::uuid const& uuid,
                    int dbLoad,
                    std::string const& dbFilename,
                    std::string const& dbId):
    super_type(name, uuid),
    M_rebuildDb(false),
    M_stage(crb::stage::online)
{
    if( ! this->findDBUuid(dbLoad, dbLoad ? dbId : dbFilename) )
        throw std::invalid_argument("Database not found during online phase");

    this->loadDB(this->absoluteDbFilename(), l );
}

template<typename RBSpace>
PBDW<RBSpace>::PBDW(std::string const& name,
                    reducedspace_ptrtype const& XR,
                    sensormap_type const& sigmas,
                    uuids::uuid const& uuid,
                    bool rebuildDb,
                    int dbLoad,
                    std::string const& dbFilename,
                    std::string const& dbId):
    super_type(name, uuid),
    M_XR(XR),
    M_sigmas(sigmas),
    M_rebuildDb(rebuildDb),
    M_stage(crb::stage::offline)
{
    this->M_M = M_sigmas.size();
    this->M_N = M_XR->size();
    for( auto const& [name,_] : M_sigmas )
        M_sensorNames.push_back(name);

    if( ! this->findDBUuid(dbLoad, dbLoad ? dbId : dbFilename) )
        this->setDBDirectory(Environment::randomUUID(true));

    if( ! M_rebuildDb )
        this->loadDB(this->absoluteDbFilename(), crb::load::all );
}

template<typename RBSpace>
typename PBDW<RBSpace>::sparse_matrix_ptrtype
PBDW<RBSpace>::initRiesz()
{
    auto u = M_XR->functionSpace()->element();
    M_qs = std::vector<element_type>(this->M_M, M_XR->functionSpace()->element());
    auto a = form2(_test=M_XR->functionSpace(), _trial=M_XR->functionSpace());
    a = integrate(_range=elements(support(M_XR->functionSpace())), _expr=inner(vf::id(u),vf::idt(u)) + inner(grad(u),gradt(u)) );
    auto am = a.matrixPtr();
    int m = 0;
    for( auto const& it /*[name, sensor]*/ : M_sigmas )
    {
        auto name = it.first;
        auto sensor = it.second;
        auto f = form1(_test=M_XR->functionSpace(), _vector=sensor->containerPtr());
        a.solve(_solution=M_qs[m++], _rhs=f, _name="pbdw");
    }
    return am;
}

template<typename RBSpace>
void
PBDW<RBSpace>::offline()
{
    if( ! M_rebuildDb )
        return;

    Feel::cout << "offline phase start with M=" << this->M_M << " and N=" << this->M_N << std::endl;
    auto u = M_XR->functionSpace()->element();
    auto am = this->initRiesz();
    this->M_matrix = matrixN_type::Zero(this->M_M+this->M_N, this->M_M+this->M_N);
    for(int i = 0; i < this->M_M; ++i )
    {
        for(int j = 0; j < i; ++j )
        {
            this->M_matrix(i, j) = am->energy(M_qs[i], M_qs[j]);
            this->M_matrix(j, i) = this->M_matrix(i, j);
        }
        this->M_matrix(i, i) = am->energy(M_qs[i], M_qs[i]);
        for(int j = 0; j < this->M_N; ++j )
            this->M_matrix(i, this->M_M+j) = am->energy(M_qs[i], M_XR->primalBasisElement(j));
    }
    this->M_matrix.bottomLeftCorner(this->M_N, this->M_M) = this->M_matrix.topRightCorner(this->M_M, this->M_N).transpose();

    Feel::cout << "computing " << this->M_Nl << " outputs" << std::endl;
    this->M_F = matrixN_type::Zero(this->M_Nl, this->M_M+this->M_N);
    for( int i = 0; i < this->M_Nl; ++i )
    {
        for( int j = 0; j < this->M_M; ++j )
            this->M_F(i,j) = inner_product( *this->M_Fs[i], this->M_qs[j] );
        for( int j = this->M_M; j < this->M_M+this->M_N; ++j )
            this->M_F(i,j) = inner_product( *this->M_Fs[i], this->M_XR->primalBasisElement(j-M_M) );
    }

    this->saveDB();
}

template<typename RBSpace>
void
PBDW<RBSpace>::setOutputs( std::vector<vector_ptrtype> const& Fs, std::vector<std::string> const& names )
{
    this->M_Fs = Fs;
    this->M_Nl = Fs.size();
    if( names.size() == this->M_Nl )
        this->M_outputNames = names;
    else
    {
        for(int i = 0; i < this->M_Nl; ++i)
            M_outputNames.push_back("output_"+std::to_string(i));
    }
}

template<typename RBSpace>
typename PBDW<RBSpace>::element_type
PBDW<RBSpace>::solution(vectorN_type const& yobs) const
{
    int M = this->dimensionM();
    int N = this->dimensionN();
    auto vn = this->online(yobs);
    // Feel::cout << "basis coeffients:\n" << vn << std::endl;
    auto wn = M_XR->primalRB();
    auto Xh = M_XR->functionSpace();
    auto I = Xh->element();
    I.zero();
    for(int i = 0; i < M; ++i )
        I.add(vn(i), M_qs[i]);
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
        j["markers"] = M_markers;
        std::ofstream o(this->absoluteJsonFilename());
        o << j.dump(2) << std::endl;

        fs::ofstream ofs( this->absoluteDbFilename() );
        if ( ofs )
        {
            Feel::cout << "saving DB at " << this->absoluteDbFilename() << std::endl;
            boost::archive::binary_oarchive oa( ofs );
            oa << this->dimensionN();
            oa << this->dimensionM();
            oa << this->M_matrix;
            oa << this->M_Nl;
            oa << this->M_F;
            oa << this->M_sensorNames;
            oa << this->M_outputNames;
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
        for( int n = 0; n < this->M_N; ++n )
            oa << M_XR->primalBasisElement(n);
        oa << M_sigmas;
        for( int m = 0; m < this->M_M; ++m )
            oa << M_qs[m];
    }
}

template<typename RBSpace>
void
PBDW<RBSpace>::loadDB( std::string const& filename, crb::load l )
{
    if( !fs::exists(filename) )
    {
        this->M_rebuildDb = true;
        return;
    }

    super_type::loadDB( filename, l);

    if( l > crb::load::rb )
    {
        std::ifstream i(this->absoluteJsonFilename());
        json j = json::parse(i);
        auto meshfilename = j["mesh"].get<std::string>();
        M_markers = j["markers"].get<std::set<std::string>>();
        auto mesh = loadMesh(_mesh=new mesh_type, _filename=meshfilename);
        space_ptrtype Xh;
        if( M_markers.empty() )
            Xh = space_type::New(_mesh=mesh, _range=elements(mesh));
        else
            Xh = space_type::New(_mesh=mesh, _range=markedelements(mesh, M_markers));
        M_XR = std::make_shared<reducedspace_type>(Xh);
        fs::ifstream ifsp( this->absoluteDbFilenameProc() );
        if( ifsp )
        {
            boost::archive::binary_iarchive ia( ifsp );
            auto u = Xh->element();
            for( int n = 0; n < this->M_N; ++n )
            {
                ia >> u;
                M_XR->addPrimalBasisElement(u);
            }
            M_sigmas = sensormap_type(Xh);
            ia >> M_sigmas;
            M_qs.resize(this->M_M);
            for( int m = 0; m < this->M_M; ++m )
            {
                ia >> u;
                M_qs[m] = Xh->element();
                M_qs[m] = u;
            }
        }
    }
}

} // namespace Feel

#endif
