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
#ifndef _FEELPP_GEIM_HPP
#define _FEELPP_GEIM_HPP 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelmor/crbdb.hpp>
#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelcore/unwrapptr.hpp>
#include <feel/feeldiscr/fsfunctionallinear.hpp>
#include <feel/feelcore/json.hpp>

#include <limits>

namespace Feel
{

/**
 * Class to support Generalized Empirical Interpolation Method.
 */
template<typename FunctionSpace>
class GEIM : public CRBDB
{
public:
    using super_type = CRBDB;
    using functionspace_type = FunctionSpace;
    using functionspace_ptrtype = std::shared_ptr<functionspace_type>;
    using mesh_type = typename functionspace_type::mesh_type;
    using element_type = typename functionspace_type::element_type;
    using linearform_type = FsFunctionalLinear<functionspace_type>;
    using linearform_ptrtype = std::shared_ptr<linearform_type>;

    using self_type = GEIM<functionspace_type>;
    using self_ptrtype = std::shared_ptr<self_type>;

    using parameterspace_type = ParameterSpaceX;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameter_type = typename parameterspace_type::element_type;
    using sampling_type = typename parameterspace_type::sampling_type;
    using sampling_ptrtype = typename parameterspace_type::sampling_ptrtype;

    using backend_type = Backend<double>;
    using vector_type = typename backend_type::vector_type;
    using vector_ptrtype = typename backend_type::vector_ptrtype;

    using solver_type = std::function<element_type (parameter_type const& mu)>;

    using vectorN_type = Eigen::VectorXd;
    using matrixN_type = Eigen::MatrixXd;

    using index_type = size_type;

public:
    /**
     * Constructor for online phase.
     * @param name Name of geim
     * @param l Loading type (rb,fe,all)
     * @param uuid Uuid to use for db
     * @param dbLoad Loading type for the DB
     * @param dbFilename Filename of the DB for load type filename
     * @param dbIf If of the DB for load type id
     */
    explicit GEIM(std::string const& name, crb::load l = crb::load::rb,
                  uuids::uuid const& uuid = uuids::nil_uuid(),
                  int dbLoad = ioption("geim.db.load"),
                  std::string const& dbFilename = soption("geim.db.filename"),
                  std::string const& dbId = soption("geim.db.id"));
    /**
     * Constructor for offline phase.
     * @param name Name of geim
     * @param sigma Linear forms
     * @param sampling Training set
     * @param solver Function to use to get an element from a parameter
     * @param uuid Uuid to use for db
     * @param rebuildDB Boolean to rebuild the database
     * @param dbLoad Loading type for the DB
     * @param dbFilename Filename of the DB for load type filename
     * @param dbIf If of the DB for load type id
     */
    GEIM(std::string const& name,
         std::vector<linearform_ptrtype> const& sigma,
         sampling_ptrtype sampling,
         solver_type solver,
         uuids::uuid const& uuid = uuids::nil_uuid(),
         bool rebuildDB = boption("geim.rebuild-database"),
         int dbLoad = ioption("geim.db.load"),
         std::string const& dbFilename = soption("geim.db.filename"),
         std::string const& dbId = soption("geim.db.id"));
    int dimension() const { return M_M; } /**< dimension of the geim */
    /**
     * Beta coefficients.
     * @param mu Parameter to use
     * @return Beta coefficients for interpolation
     */
    vectorN_type beta( parameter_type const& mu );
    /**
     * Beta coefficients.
     * @param vn Values of linear forms for interpolation
     * @return Beta coefficients for interpolation
     */
    vectorN_type beta( vectorN_type const& vn );
    /**
     * Beta coefficients.
     * @param v Element to use
     * @return Beta coefficients for interpolation
     */
    vectorN_type beta( element_type const& v);
    /**
     * Interpolation of solver(mu)
     * @param mu Parameter to use
     * @return Interpolation
     */
    element_type interpolant( parameter_type const& mu );
    /**
     * Interpolation for u
     * @param u Element to use
     * @return Interpolation
     */
    element_type interpolant( element_type const& u );
    /**
     * Interpolation for element corresponding to evaluations
     * @param vn Evaluations of linear forms
     * @return Beta coefficients for interpolation
     */
    element_type interpolant( vectorN_type const& vn );
    matrixN_type matrixB() const { return M_B; } /**< Interpolation matrix */
    std::vector<element_type> q() const { return M_q; } /**< Interpolation basis */
    element_type q( int i ) const { return M_q[i]; } /**< i-th interpolation basis */
    std::vector<parameter_type> mus() const { return M_mus; } /**< Parameters used */
    std::vector<linearform_ptrtype> sigmas() const { return M_sigmas; } /**< Set of linear forms */
    std::vector<index_type> indices() const { return M_indices; } /**< Set of indices of linear forms used */
    functionspace_ptrtype space() const { return M_Xh; }
    void setSolver( solver_type f ) { M_solver = f; } /**< Set solver to use */
    void setMaxM( int M ) { M_MaxM = M; } /**< Set maximum number of basis */
    void setTolerance( double tol ) { M_tol = tol; } /**< Set tolerance */
    void offline(); /**< Do offline phase */

protected:
    virtual index_type maxError();
    virtual double applyForm(index_type i, element_type const& u);
    virtual element_type getElement(parameter_type const& mu) { return M_solver(mu); }
    void loadDB( std::string const& filename, crb::load l ) override;
    void saveDB() override;

private:
    std::string M_name;
    functionspace_ptrtype M_Xh;
    std::vector<linearform_ptrtype> M_sigmas;
    std::vector<index_type> M_indices;
    sampling_ptrtype M_sampling;
    std::vector<parameter_type> M_mus;
    std::vector<element_type> M_q;
    matrixN_type M_B;
    int M_MaxM;
    int M_M;
    solver_type M_solver;
    double M_tol;

    element_type M_currentQ;

    bool M_rebuildDb;
    int M_dbLoad;
    std::string M_dbFilename;
    std::string M_dbId;
    crb::stage M_stage;
};

template<typename FunctionSpace>
GEIM<FunctionSpace>::GEIM(std::string const& name, crb::load l,
                          uuids::uuid const& uuid,
                          int dbLoad,
                          std::string const& dbFilename,
                          std::string const& dbId):
    super_type(name, "geim", uuid),
    M_name(name),
    M_rebuildDb(false),
    M_dbLoad(dbLoad),
    M_dbFilename(dbFilename),
    M_dbId(dbId),
    M_stage(crb::stage::online)
{
    if( ! this->findDBUuid(M_dbLoad, M_dbLoad ? M_dbId : M_dbFilename) )
        throw std::invalid_argument("Database not found during online phase");

    this->loadDB(this->absoluteDbFilename(), l );
}

template<typename FunctionSpace>
GEIM<FunctionSpace>::GEIM(std::string const& name,
                          std::vector<linearform_ptrtype> const& sigma,
                          sampling_ptrtype sampling,
                          solver_type solver,
                          uuids::uuid const& uuid,
                          bool rebuildDb,
                          int dbLoad,
                          std::string const& dbFilename,
                          std::string const& dbId):
    super_type( name, "geim", uuid ),
    M_name(name),
    M_sigmas(sigma),
    M_sampling(sampling),
    M_solver(solver),
    M_MaxM(ioption("geim.dimension-max")),
    M_M(0),
    M_tol(doption("geim.tolerance")),
    M_rebuildDb(rebuildDb),
    M_dbLoad(dbLoad),
    M_dbFilename(dbFilename),
    M_dbId(dbId),
    M_stage(crb::stage::offline)
{
    if( M_sigmas.size() == 0 )
        throw std::invalid_argument("No linear forms provided");

    if( ! this->findDBUuid(M_dbLoad, M_dbLoad ? M_dbId : M_dbFilename) )
        this->setDBDirectory(Environment::randomUUID(true));

    M_currentQ = M_solver(get<0>(M_sampling->max()));
    M_Xh = M_currentQ.functionSpace();

    if( M_rebuildDb )
    {
        M_M = 0;
        M_B = matrixN_type(0,0);
    }
    else
        this->loadDB(this->absoluteDbFilename(), crb::load::all );
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::vectorN_type
GEIM<FunctionSpace>::beta( parameter_type const& mu )
{
    vectorN_type rhs(M_M);
    if( M_M == 0 )
        return rhs;

    auto phi = this->getElement(mu);
    return this->beta(phi);
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::vectorN_type
GEIM<FunctionSpace>::beta( vectorN_type const& vn )
{
    vectorN_type b(M_M);
    if( M_M == 0 )
        return b;

    b = M_B.colPivHouseholderQr().solve(vn);
    return b;
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::vectorN_type
GEIM<FunctionSpace>::beta( element_type const& v )
{
    vectorN_type rhs(M_M);
    if( M_M == 0 )
        return rhs;

    for(int i = 0; i < M_M; ++i )
        rhs(i) = this->applyForm(M_indices[i], v);
    vectorN_type b = M_B.colPivHouseholderQr().solve(rhs);
    return b;
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::element_type
GEIM<FunctionSpace>::interpolant(parameter_type const& mu)
{
    auto u = this->getElement(mu);
    return this->interpolant(u);
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::element_type
GEIM<FunctionSpace>::interpolant(element_type const& u)
{
    auto I = M_Xh->element();
    I.zero();
    auto b = this->beta(u);
    for(int i = 0; i < M_M; ++i)
        I.add(b(i), M_q[i]);
    return I;
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::element_type
GEIM<FunctionSpace>::interpolant(vectorN_type const& vn)
{
    auto I = M_Xh->element();
    I.zero();
    auto b = this->beta(vn);
    for(int i = 0; i < M_M; ++i)
        I.add(b(i), M_q[i]);
    return I;
}

template<typename FunctionSpace>
typename GEIM<FunctionSpace>::index_type
GEIM<FunctionSpace>::maxError()
{
    int maxI = 0;
    double max = std::numeric_limits<double>::lowest();
    for( int i = 0; i < this->M_sigmas.size(); ++i )
    {
        // skip existing forms
        if( std::count(M_indices.begin(), M_indices.end(), i) )
            continue;
        double c = std::abs(this->applyForm(i, M_currentQ));
        if( c > max )
        {
            max = c;
            maxI = i;
        }
    }
    return maxI;
}

template<typename FunctionSpace>
double
GEIM<FunctionSpace>::applyForm(index_type i, element_type const& u)
{
    return inner_product(M_sigmas[i]->container(), u);
}

template<typename FunctionSpace>
void
GEIM<FunctionSpace>::offline()
{
    if(!M_solver)
    {
        Feel::cout << "You need to set a solver before the offline phase!" << std::endl;
        return;
    }

    if( M_sigmas.size() > 0 )
    {
        if( M_MaxM > M_sampling->size() || M_MaxM > M_sigmas.size() )
            M_MaxM = M_sampling->size() > M_sigmas.size() ? M_sigmas.size() : M_sampling->size();
    }
    else
    {
        if( M_MaxM > M_sampling->size() )
            M_MaxM = M_sampling->size();
    }

    parameter_type mu;
    double error;

    Feel::cout << "offline phase starts" << std::endl;
    while( M_M < M_MaxM )
    {
        double nMax = std::numeric_limits<double>::lowest();
        for( auto const& m : *M_sampling )
        {
            auto phi = this->getElement(m);
            auto I = this->interpolant(phi);
            phi -= I;
            auto n = phi.l2Norm();
            if( n > nMax )
            {
                nMax = n;
                mu = m;
            }
        }

        M_currentQ = this->getElement(mu);
        auto I = this->interpolant(M_currentQ);
        M_currentQ -= I;
        error = M_currentQ.l2Norm();
        Feel::cout << "error with " << M_M << " basis = " << error << std::endl;
        if( error < M_tol )
            break;
        auto index = this->maxError();
        M_currentQ.scale(1./this->applyForm(index, M_currentQ));

        M_mus.push_back(mu);
        M_indices.push_back(index);
        M_q.push_back(M_currentQ);

        M_B.conservativeResize(M_M+1, M_M+1);
        for( int i = 0; i < M_M; ++i )
        {
            M_B(M_M, i) = this->applyForm(index, M_q[i]);
            M_B(i, M_M) = this->applyForm(M_indices[i], M_currentQ); // should be 0
        }
        M_B(M_M, M_M) = 1.;
        M_M++;

        this->saveDB();
    }
}

template<typename FunctionSpace>
void
GEIM<FunctionSpace>::saveDB()
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
        o << j.dump(2) << std::endl;
        fs::ofstream ofs( this->absoluteDbFilename() );
        if ( ofs )
        {
            Feel::cout << "saving DB at " << this->absoluteDbFilename() << std::endl;
            //boost::archive::text_oarchive oa( ofs );
            boost::archive::binary_oarchive oa( ofs );
            // write class instance to archive
            oa << M_M;
            oa << M_indices;
            oa << M_B;
            oa << M_mus;
        }
    }
    if( ! fs::exists(this->absoluteMeshFilename()) )
    {
        Feel::cout << "Saving mesh to " << this->absoluteMeshFilename() << std::endl;
        M_Xh->mesh()->saveHDF5( this->absoluteMeshFilename() );
    }
    fs::ofstream ofsp( this->absoluteDbFilenameProc() );
    if( ofsp )
    {
        boost::archive::binary_oarchive oa( ofsp );
        oa << M_q;
        // oa << M_sigmas.size();
        // oa << M_sampling;
        // oa << M_sigmas;
        // archive and stream closed when destructors are called
    }
}

template<typename FunctionSpace>
void
GEIM<FunctionSpace>::loadDB( std::string const& filename, crb::load l )
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
            ia >> M_M;
            M_indices.resize(M_M);
            ia >> M_indices;
            M_B.resize(M_M,M_M);
            ia >> M_B;
        }
        if( l > crb::load::fe )
        {
            M_mus.resize(M_M);
            ia >> M_mus;
            // int sigmaSize;
            // ia >> sigmaSize;
            // M_sigmas.resize(sigmaSize);
            // ia >> M_sampling;
            // ia >> M_sigmas;
        }
        // archive and stream closed when destructors are called
    }
    if( l > crb::load::rb )
    {
        if( ! M_Xh )
        {
            std::ifstream i(this->absoluteJsonFilename());
            json j = json::parse(i);
            auto meshfilename = j["mesh"].get<std::string>();
            auto mesh = loadMesh(_mesh=new mesh_type, _filename=meshfilename);
            M_Xh = functionspace_type::New(mesh);
        }
        fs::ifstream ifsp( this->absoluteDbFilenameProc() );
        if( ifsp )
        {
            boost::archive::binary_iarchive ia( ifsp );
            M_q.resize(M_M);
            for(int i = 0; i < M_M; ++i )
                M_q[i] = M_Xh->element();
            ia >> M_q;
        }
    }
}

} // namespace Feel

#endif
