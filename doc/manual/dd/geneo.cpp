/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
             Pierre Jolivet <pierre.jolivet@imag.fr>
       Date: 2013-12-21

  Copyright (C) 2013 Université de Strasbourg
  Copyright (C) 2013 Université de Grenoble

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
/**
   \file geneo.cpp
   \author Pierre Jolivet <pierre.jolivet@imag.fr>
   \date 2013-12-21
 */
#ifdef EIGEN_USE_MKL_ALL
#undef EIGEN_USE_MKL_ALL
#endif
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelalg/vectoreigen.hpp>
#include <feel/feelalg/matrixeigensparse.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feeltiming/tic.hpp>
#ifdef FEELPP_HAS_HPDDM
#define BDD            // BDD module
// #define FETI           // FETI module
#define MUMPSSUB       // MUMPS as solver inside subdomain
#define DMUMPS         // MUMPS as distributed solver for the coarse operator
#include <HPDDM.hpp>
#endif

#include "geneo_coefficients.hpp"

namespace Feel
{
template<uint16_type Dim, uint16_type Order>
static inline void generateRBM(unsigned short& nb, double**& ev, boost::shared_ptr<FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Scalar>>>>& Vh) {
    nb = 1;
    ev = new double*[nb];
    *ev = new double[nb * Vh->nDof()];
    auto rbm = Vh->element();
    rbm = vf::project(Vh, elements(Vh->mesh()), cst(1.0));
    std::copy(rbm.begin(), rbm.end(), *ev);
}
template<uint16_type Order>
static inline void generateRBM(unsigned short& nb, double**& ev, boost::shared_ptr<FunctionSpace<Mesh<Simplex<2>>, bases<Lagrange<Order, Vectorial>>>>& Vh) {
    nb = 3;
    ev = new double*[nb];
    *ev = new double[nb * Vh->nDof()];
    for(unsigned short i = 0; i < nb; ++i)
        ev[i] = *ev + i * Vh->nDof();
    auto rbm = Vh->element();
    rbm = vf::project(Vh, elements(Vh->mesh()), oneX());
    std::copy(rbm.begin(), rbm.end(), ev[0]);
    rbm = vf::project(Vh, elements(Vh->mesh()), oneY());
    std::copy(rbm.begin(), rbm.end(), ev[1]);
    rbm = vf::project(Vh, elements(Vh->mesh()), vec(Py(), -Px()));
    std::copy(rbm.begin(), rbm.end(), ev[2]);
}
template<uint16_type Order>
static inline void generateRBM(unsigned short& nb, double**& ev, boost::shared_ptr<FunctionSpace<Mesh<Simplex<3>>, bases<Lagrange<Order, Vectorial>>>>& Vh) {
    nb = 6;
    ev = new double*[nb];
    *ev = new double[nb * Vh->nDof()];
    for(unsigned short i = 0; i < nb; ++i)
        ev[i] = *ev + i * Vh->nDof();
    auto rbm = Vh->element();
    rbm = vf::project(Vh, elements(Vh->mesh()), oneX());
    std::copy(rbm.begin(), rbm.end(), ev[0]);
    rbm = vf::project(Vh, elements(Vh->mesh()), oneY());
    std::copy(rbm.begin(), rbm.end(), ev[1]);
    rbm = vf::project(Vh, elements(Vh->mesh()), oneZ());
    std::copy(rbm.begin(), rbm.end(), ev[2]);
    rbm = vf::project(Vh, elements(Vh->mesh()), vec(Py(), -Px(), cst(0.0)));
    std::copy(rbm.begin(), rbm.end(), ev[3]);
    rbm = vf::project(Vh, elements(Vh->mesh()), vec(-Pz(), cst(0.0), Px()));
    std::copy(rbm.begin(), rbm.end(), ev[4]);
    rbm = vf::project(Vh, elements(Vh->mesh()), vec(cst(0.0), Pz(), -Py()));
    std::copy(rbm.begin(), rbm.end(), ev[5]);
}
template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type>
static inline void coefficients(double* r, boost::shared_ptr<FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>>>>& Vh) {
    auto x = Vh->element();
    if(!std::is_same<Type<Dim>, Vectorial<Dim>>::value) {
        kappa k;
        k.val = doption("parameters.kappa");
        x = vf::project(Vh, elements(Vh->mesh()), idf(k) * one());
    }
    else {
        stripes s;
        s.first  = doption("parameters.n");
        s.second = doption("parameters.nu");
        x = vf::project(Vh, elements(Vh->mesh()), idf(s) * one());
    }
    std::copy(x.begin(), x.end(), r);
}
template<uint16_type Dim, uint16_type Order>
static inline void assemble(Backend<double>::sparse_matrix_ptrtype& A, Backend<double>::vector_ptrtype& f, boost::shared_ptr<FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Scalar>>>>& Vh) {
    typename FunctionSpace<typename Mesh<Simplex<Dim>>::mesh_type, bases<Lagrange<Order, Scalar>>>::element_type u = Vh->element(), v = Vh->element();
    kappa k;
    k.val = doption("parameters.kappa");
    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(Vh->mesh()), _expr = idf(k) * gradt(u) * trans(grad(v)));
    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(Vh->mesh()), _expr = id(v));
    a += on(_range = markedfaces(Vh->mesh(), "Dirichlet"), _rhs = l, _element = u, _expr = cst(0.0));
    A->close();
}
template<uint16_type Order>
static inline void assemble(Backend<double>::sparse_matrix_ptrtype& A, Backend<double>::vector_ptrtype& f, boost::shared_ptr<FunctionSpace<Mesh<Simplex<2>>, bases<Lagrange<Order, Vectorial>>>>& Vh) {
    typename FunctionSpace<typename Mesh<Simplex<2>>::mesh_type, bases<Lagrange<Order, Vectorial>>>::element_type u = Vh->element(), v = Vh->element();
    stripes s;
    s.first  = doption("parameters.epsilon");
    s.second = doption("parameters.e");
    auto E = idf(s);
    s.first  = doption("parameters.n");
    s.second = doption("parameters.nu");
    auto nu = idf(s);
    auto mu = E / (2 * (1 + nu));
    auto lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

    const double density = 1.0e+3;

    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(Vh->mesh()),
                  _expr = lambda * divt(u) * div(v) +
                          2 * mu * trace(trans(sym(gradt(u))) * sym(grad(u))));
    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(Vh->mesh()), _expr = -density * trans(oneY()) * id(v));
    a += on(_range = markedfaces(Vh->mesh(), "Dirichlet"), _rhs = l, _element = u, _expr = zero<2, 1>());
    A->close();
}
template<uint16_type Order>
static inline void assemble(Backend<double>::sparse_matrix_ptrtype& A, Backend<double>::vector_ptrtype& f, boost::shared_ptr<FunctionSpace<Mesh<Simplex<3>>, bases<Lagrange<Order, Vectorial>>>>& Vh) {
    typename FunctionSpace<typename Mesh<Simplex<3>>::mesh_type, bases<Lagrange<Order, Vectorial>>>::element_type u = Vh->element(), v = Vh->element();
    stripes s;
    s.first  = doption("parameters.epsilon");
    s.second = doption("parameters.e");
    auto E = idf(s);
    s.first  = doption("parameters.n");
    s.second = doption("parameters.nu");
    auto nu = idf(s);
    auto mu = E / (2 * (1 + nu));
    auto lambda = E * nu / ((1 + nu) * (1 - 2 * nu));

    const double density = 1.0e+3;

    auto a = form2(_trial = Vh, _test = Vh, _matrix = A);
    a = integrate(_range = elements(Vh->mesh()),
                  _expr = lambda * divt(u) * div(v) +
                          2 * mu * trace(trans(sym(gradt(u))) * sym(grad(u))));
    auto l = form1(_test = Vh, _vector = f);
    l = integrate(_range = elements(Vh->mesh()), _expr = -density * trans(vec(cst(1.0/std::sqrt(2.0)), cst(0.0), cst(1.0/std::sqrt(2.0)))) * id(v));
    a += on(_range = markedfaces(Vh->mesh(), "Dirichlet"), _rhs = l, _element = u, _expr = zero<3, 1>());
    A->close();
}

template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type>
class Geneopp : public Simget
{
public:
    void run();
}; // Geneopp

template<uint16_type Dim, uint16_type Order, template<uint16_type> class Type>
void Geneopp<Dim, Order, Type>::run()
{
    static_assert(Dim == 2 || Dim == 3, "Wrong dimension");
#if defined(FETI)
    typedef HpFeti<FetiPrcdtnr::DIRICHLET, double, 'S'> prec_type;
    constexpr unsigned short nu = 0;
#elif defined(BDD)
    typedef HpBdd<double, 'S'> prec_type;
    unsigned short nu = ioption("nu");
#endif
    prec_type* K = new prec_type;
    tic();
    unsigned short p = ioption("p");
    unsigned short topology = ioption("topology");
    MPI_Comm comm;
    bool exclude = boption("exclude");
    bool excluded = HPDDM::DMatrix::splitCommunicator(Environment::worldComm(), comm, exclude, p, topology);
    if(!exclude)
        MPI_Comm_dup(Environment::worldComm(), &comm);
    boost::mpi::communicator bComm(comm, boost::mpi::comm_take_ownership);
    WorldComm wComm(bComm, bComm, bComm, bComm.rank(), std::vector<int>(bComm.size(), true));
    toc("communicators");
    std::vector<double> timers(6 + (nu != 0 ? 2 : 0) + (exclude == true));
    mpi::timer time;
    boost::shared_ptr<Mesh<Simplex<Dim>>> mesh;
    typename FunctionSpace<typename Mesh<Simplex<Dim>>::mesh_type, bases<Lagrange<Order, Type>>>::element_type uLocal;
    double* b;
    std::vector<int> interface;
    if(!excluded) {
        tic();
        if(std::is_same<Type<Dim>, Vectorial<Dim>>::value) {
            mesh = createGMSHMesh(_mesh = new Mesh<Simplex<Dim>>(wComm),
                                  _update = MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                  _worldcomm = wComm,
                                  _desc = domain(_worldcomm = wComm, _name = "hypercube", _shape = "hypercube",
                                                 _xmin = 0.0, _xmax = 10.0,
                                                 _ymin = 0.0, _ymax = 1.0,
                                                 _zmin = 0.0, _zmax = 1.0));
            mesh->addMarkerName("Dirichlet", Dim == 2 ? 1 : 19, Dim == 2 ? 1 : 2);
        }
        else {
            mesh = createGMSHMesh(_mesh = new Mesh<Simplex<Dim>>(wComm),
                                  _update = MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                  _worldcomm = wComm,
                                  _desc = domain(_worldcomm = wComm, _name = "hypercube", _shape = "hypercube",
                                                 _xmin = 0.0, _xmax = 1.0,
                                                 _ymin = 0.0, _ymax = 1.0,
                                                 _zmin = 0.0, _zmax = 1.0));
            mesh->addMarkerName("Dirichlet", Dim == 2 ? 1 : 19, Dim == 2 ? 1 : 2);
        }
        toc("global mesh");
        bComm.barrier();
        time.restart();
        double** ev;
        auto meshLocal = createSubmesh(mesh, elements(mesh), Environment::worldCommSeq());
        timers[0] = time.elapsed();
        bComm.barrier();
        time.restart();
        auto VhLocal = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>>>::New(_mesh = meshLocal,
                       _worldscomm = Environment::worldsCommSeq(1));
        timers[1] = time.elapsed();
        bComm.barrier();
        time.restart();
        {
            std::vector<std::vector<int>> map(mesh->faceNeighborSubdomains().size());
            std::string kind = soption(_name = "backend");
            boost::shared_ptr<Backend<double>> ptr_backend = Backend<double>::build(kind, "", Environment::worldCommSeq());
            Backend<double>::vector_ptrtype f = ptr_backend->newVector(VhLocal);
#pragma omp parallel for shared(map) schedule(static, 4)
            for(int i = 0; i < map.size(); ++i) {
                auto trace = createSubmesh(mesh, interprocessfaces(mesh, *std::next(mesh->faceNeighborSubdomains().begin(), i)),
                                           Environment::worldCommSeq());
                auto Xh = FunctionSpace<typename Mesh<Simplex<Dim>>::trace_mesh_type, bases<Lagrange<Order, Type>>>::New(_mesh = trace, _worldscomm = Environment::worldsCommSeq(1));
                auto l = ptr_backend->newVector(Xh);
                double* pt = &(*l)(0);
                std::iota(pt, pt + Xh->nDof(), 1.0);
                auto op = opInterpolation(_domainSpace = VhLocal,
                                          _imageSpace  = Xh,
                                          _backend     = ptr_backend, _ddmethod = true);
                map[i].resize(Xh->nDof());
#pragma omp critical
                {
                    ptr_backend->prod(op->mat(), *l, *f, true);
                    pt = &(*f)(0);
                    for(int j = 0; j < VhLocal->nDof(); ++j, ++pt)
                        if(std::round(*pt) != 0)
                            map[i][std::round(*pt) - 1] = j;
                }
            }
            {
                std::set<int> unique;
                for(std::vector<int>& pt : map)
                    for(int& i : pt)
                        unique.insert(i);
                interface.insert(interface.begin(), unique.cbegin(), unique.cend());
                std::unordered_map<int, int> mapping;
                mapping.reserve(interface.size());
                int j = 0;
                for(const int& i : interface)
                    mapping[i] = j++;
                for(std::vector<int>& pt : map)
                    for(int& i : pt)
                        i = mapping[i];
            }
            timers[2] = time.elapsed();
            bComm.barrier();
            time.restart();
            Backend<double>::sparse_matrix_ptrtype A = ptr_backend->newMatrix(VhLocal, VhLocal);
            assemble(A, f, VhLocal);
            timers[3] = time.elapsed();
            b = new double[VhLocal->nDof()];
#ifdef FEELPP_HAS_HPDDM
            int* ic = new int[VhLocal->nDof() + 1];
            double* c;
            int* jc;
            if(kind == "petsc") {
                Mat PetscA = static_cast<MatrixPetsc<double>*>(&*A)->mat();
                Vec PetscF = static_cast<VectorPetsc<double>*>(&*f)->vec();
                PetscInt n;
                const PetscInt* ia;
                const PetscInt* ja;
                PetscScalar* array;
                PetscBool done;
                MatGetRowIJ(PetscA, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
                MatSeqAIJGetArray(PetscA, &array);
                std::copy_n(ia, n + 1, ic);
                c = new double[ia[n]];
                jc = new int[ia[n]];
                std::copy_n(array, ia[n], c);
                std::copy_n(ja, ia[n], jc);
                MatSeqAIJRestoreArray(PetscA, &array);
                MatRestoreRowIJ(PetscA, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
                VecGetArray(PetscF, &array);
                std::copy_n(array, VhLocal->nDof(), b);
                VecRestoreArray(PetscF, &array);
            }
            else {
                Eigen::SparseMatrix<double>& EigenA = static_cast<MatrixEigenSparse<double>*>(&*A)->mat();
                EigenA.makeCompressed();
                Eigen::Matrix<double, Eigen::Dynamic, 1>& EigenF = static_cast<VectorEigen<double>*>(&*f)->vec();
#if 0
                Eigen::SparseMatrix<double, Eigen::RowMajor> EigenC(EigenA);
#endif
                std::copy_n(EigenA.outerIndexPtr(), VhLocal->nDof() + 1, ic);
                c = new double[ic[VhLocal->nDof()]];
                jc = new int[ic[VhLocal->nDof()]];
                std::copy_n(EigenA.innerIndexPtr(), ic[VhLocal->nDof()], jc);
                std::copy_n(EigenA.valuePtr(), ic[VhLocal->nDof()], c);
                std::copy_n(EigenF.data(), EigenF.rows(), b);
            }
            K->Subdomain::initialize(new HPDDM::MatrixCSR<double>(VhLocal->nDof(), VhLocal->nDof(), ic[VhLocal->nDof()], c, ic, jc, false, true), mesh->faceNeighborSubdomains(), map, &comm);
            decltype(map)().swap(map);
            if(nu == 0) {
                if(nelements(markedfaces(meshLocal, "Dirichlet")) == 0) {
                    unsigned short nb;
                    generateRBM(nb, ev, VhLocal);
                    K->setVectors(ev);
                    K->super::super::initialize(nb);
                }
                else
                    K->super::super::initialize(0);
            }
#endif
        }
#ifdef FEELPP_HAS_HPDDM
        bComm.barrier();
        time.restart();
        K->renumber(interface, b);
        timers[4] = time.elapsed();
        if(nu != 0) {
            bComm.barrier();
            time.restart();
            K->computeSchurComplement();
            timers[6] = time.elapsed();
            bComm.barrier();
            time.restart();
            K->solveGEVP<'S'>(nu);
            timers[7] = time.elapsed();
            K->super::super::initialize(nu);
        }
        K->callNumfactPreconditioner();
        const char scaling = soption("scaling")[0];
        double* rho = nullptr;
        if(scaling == 'r') {
            double* rho = new double[VhLocal->nDof()];
            coefficients(rho, VhLocal);
            K->renumber(interface, rho);
        }
        K->buildScaling(scaling, rho);
        delete [] rho;
        uLocal = VhLocal->element();
    }
    std::vector<unsigned short> parm(5);
    parm[HPDDM::Parameter::P]            = p;
    parm[HPDDM::Parameter::TOPOLOGY]     = topology;
    parm[HPDDM::Parameter::DISTRIBUTION] = HPDDM::DMatrix::NON_DISTRIBUTED;
    parm[HPDDM::Parameter::STRATEGY]     = ioption("strategy");
    if(nu == 0)
        parm[HPDDM::Parameter::NU]       = excluded ? 0 : K->getLocal();
    else
        parm[HPDDM::Parameter::NU]       = nu;
    unsigned short it = ioption("it");
    double eps = doption("eps");
    if(excluded) {
        K->Subdomain::initialize(&comm);
        K->buildTwo<2>(Environment::worldComm(), parm);
        Environment::worldComm().barrier();
        HPDDM::IterativeMethod::PCG<true>(*K, static_cast<double*>(nullptr), static_cast<double*>(nullptr), it, eps, Environment::worldComm(), Environment::isMasterRank());
        delete K;
    }
    else {
        std::pair<MPI_Request, const double*>* ret;
        int flag = 0;
        if(exclude) {
            ret = K->buildTwo<1>(Environment::worldComm(), parm);
            if(ret) {
                MPI_Test(&(ret->first), &flag, MPI_STATUS_IGNORE);
                if(flag) {
                    delete [] ret->second;
                    delete ret;
                }
            }
        }
        else
            ret = K->buildTwo<0>(Environment::worldComm(), parm);
        time.restart();
        K->callNumfact();
        timers[5] = time.elapsed();
        if(ret) {
            time.restart();
            if(flag)
                timers.back() = 0.0;
            else {
                MPI_Wait(&(ret->first), MPI_STATUS_IGNORE);
                delete [] ret->second;
                delete ret;
            }
            timers.back() = time.elapsed();
            time.restart();
        }
        Environment::worldComm().barrier();
        if(ret)
            timers.back() += time.elapsed();
        HPDDM::IterativeMethod::PCG<false>(*K, &(uLocal[0]), b, it, eps, Environment::worldComm(), Environment::isMasterRank());

        double* storage = new double[2];
        K->computeError(&(uLocal[0]), b, storage);
        delete [] b;

        K->originalNumbering(interface, &(uLocal[0]));

        double stats[3] = { K->getMult() / 2.0, mesh->faceNeighborSubdomains().size() / static_cast<double>(bComm.size()), static_cast<double>(K->getAllDof()) };
        if(bComm.rank() == 0) {
            std::streamsize ss = std::cout.precision();
            std::cout << std::scientific << " --- error = " << storage[1] << " / " << storage[0] << std::endl;
            MPI_Reduce(MPI_IN_PLACE, stats, 3, MPI_DOUBLE, MPI_SUM, 0, comm);
            std::cout << std::fixed << " --- number of Lagrange multipliers: " << static_cast<unsigned long>(stats[0]) << ", number of dof: " << static_cast<unsigned long>(stats[2]) << ", on average, number of neighbors: " << std::setprecision(1) << stats[1] << " and h: " << std::scientific << mesh->hAverage() << std::endl;
            std::cout << std::fixed;
            std::cout.precision(ss);
        }
        else
            MPI_Reduce(stats, NULL, 3, MPI_DOUBLE, MPI_SUM, 0, comm);
        delete [] storage;
#else
        uLocal = vf::project(VhLocal, elements(meshLocal), Px() + Py());
#endif
        delete K;
        if(boption("export")) {
            tic();
            time.restart();
            auto VhVisu = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order, Type>>>::New(_mesh = mesh, _worldscomm = worldsComm(wComm));
            double timeFeel = time.elapsed();
            auto uVisu = vf::project(_space = VhVisu, _expr = idv(uLocal));
            auto e = exporter(_mesh = mesh);
            e->add("u", uVisu);
            e->add("rank", regionProcess(Pdh<0>(mesh)));
            if(ioption("v") > 2) {
                time.restart();
                boost::shared_ptr<Backend<double>> ptr_global = Backend<double>::build("petsc", "", wComm);
                auto D = ptr_global->newMatrix(VhVisu, VhVisu);
                auto F = ptr_global->newVector(VhVisu);
                assemble(D, F, VhVisu);
                auto u = VhVisu->element();
                ptr_global->solve(_matrix = D, _solution = u, _rhs = F);
                timeFeel += time.elapsed();
                e->add("u_ref", u);
                u -= uVisu;
                double error = u.l2Norm() / F->l2Norm();
                if(bComm.rank() == 0)
                    std::cout << std::scientific << " --- relative error with respect to solution from PETSc = " << error << " (computed in " << timeFeel << ")" << std::endl;
            }
            e->save();
            toc("exporter");
        }
        double* allTimers = new double[timers.size() * bComm.size()];
        MPI_Gather(timers.data(), timers.size(), MPI_DOUBLE, allTimers, timers.size(), MPI_DOUBLE, 0, bComm);
        if(bComm.rank() == 0) {
            std::cout << "createSubmesh  FunctionSpace  Interface    Assem. solver  Renumbering    Fact. pinv.  ";
            if(timers.size() >= 8)
                std::cout << "Schur compl.   GenEO       ";
            if(timers.size() == 7 || timers.size() == 9)
                std::cout << "  Waiting time";
            std::cout << std::endl;
            std::cout.precision(4);
            for(int i = 0; i < bComm.size(); ++i) {
                std::cout << std::scientific << allTimers[0 + i * timers.size()] << "     " << allTimers[1 + i * timers.size()] << "     " << allTimers[2 + i * timers.size()] << "   " << allTimers[3 + i * timers.size()] << "     " << allTimers[4 + i * timers.size()] << "     " << allTimers[5 + i * timers.size()] << "   ";
                if(timers.size() >= 8)
                    std::cout << std::scientific << allTimers[6 + i * timers.size()] << "     " << allTimers[7 + i * timers.size()] << "  ";
                if(timers.size() == 7 || timers.size() == 9)
                    std::cout << "  " << std::scientific << allTimers[(i + 1) * timers.size() - 1];
                std::cout << std::endl;
            }
        }
        delete [] allTimers;
    }
} // Geneopp::run
} // Feel

int main(int argc, char** argv) {
    using namespace Feel;

    po::options_description opts("FETI/BDD options");
    opts.add_options()
        ("scaling", po::value<std::string>()->default_value("m"), "kind of scaling")
        ("p", po::value<int>()->default_value(1), "number of master processes")
        ("topology", po::value<int>()->default_value(0), "distribution of the coarse operator")
        ("strategy", po::value<int>()->default_value(3), "ordering tool for the direct coarse solver (only useful when using MUMPS)")
        ("eps", po::value<double>()->default_value(1e-8), "relative preconditioned residual")
        ("it", po::value<int>()->default_value(50), "maximum number of iterations")
        ("nu", po::value<int>()->default_value(0), "number of eigenvalues")
        ("exclude", po::value<bool>()->default_value(false), "exclude the master processes")
        ;
    /**
     * Initialize Feel++ Environment
     */
    Environment env(_argc = argc, _argv = argv, _desc = opts,
                    _about = about(_name = "geneo",
                                   _author = "Feel++ Consortium",
                                   _email = "feelpp-devel@feelpp.org"));
    Application app;
    // app.add(new Geneopp<2, 1, Scalar>());
    // app.add(new Geneopp<3, 2, Scalar>());
    app.add(new Geneopp<2, 1, Vectorial>());
    // app.add(new Geneopp<2, 7, Vectorial>());
    // app.add(new Geneopp<3, 2, Vectorial>());
    app.run();
}
