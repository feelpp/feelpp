/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file dist2wallsoptimized.cpp
   \author Guillaume Dolle <gdolle at unistra.fr>
   \date 2014-01-21
 */

//#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE precAFP
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feel.hpp>
#include <feel/feelpde/preconditionerblockms.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feeldiscr/ned1h.hpp>

#if FEELPP_DIM == 2
#define curl_op curlx
#define curlt_op curlxt
#define curlv_op curlxv
#else
#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv
#endif

using namespace Feel;


inline
po::options_description
makeOptions()
{
    po::options_description opts( "test_precAFP" );
    opts.add_options()
    ( "test-case", po::value<int>()->default_value( -1 ), "The test case number" )
    ( "title", po::value<std::string>()->default_value( "noTitle" ), "The title for jekyll" )
    ( "generateMD", po::value<bool>()->default_value( false ), "Save MD file" )
    ( "saveTimers", po::value<bool>()->default_value( true ), "print timers" )
    ( "myModel", po::value<std::string>()->default_value( "model.mod" ), "name of the model" )
    ( "penaldir", po::value<double>()->default_value( 10 ), "penalisation term value" )
    ;
    return opts.add( Feel::feel_options() )
        .add(Feel::backend_options("ms"))
        .add(Feel::blockms_options("ms"));
}

inline
AboutData
makeAbout()
{
#if FEELPP_DIM==2
    AboutData about( "precAFP2D" ,
                     "precAFP2D" ,
                     "0.1",
                     "test precAFP2D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#else
    AboutData about( "precAFP3D" ,
                     "precAFP3D" ,
                     "0.1",
                     "test precAFP3D",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );
#endif

    about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );

    return about;
}

///     \tparam DIM         Topological dimension.
template<int DIM>
class TestPrecAFP : public Application
{
    private:
    typedef Application super;
    //! Numerical type is double
    typedef double value_type;

    //! Simplexes of order ORDER
    typedef Simplex<DIM> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //! Hcurl space
    typedef Nedelec<0,NedelecKind::NED1 > curl_basis_type;
    typedef FunctionSpace<mesh_type, bases<curl_basis_type>> curl_space_type;
    typedef std::shared_ptr<curl_space_type> curl_space_ptrtype;
    typedef typename curl_space_type::element_type curl_element_type;

    //! Pch space
    typedef Lagrange<1, Scalar> lag_basis_type; 
    typedef FunctionSpace<mesh_type, bases<lag_basis_type>> lag_space_type;
    typedef std::shared_ptr<lag_space_type> lag_space_ptrtype;
    typedef typename lag_space_type::element_type lag_element_type;

    //! Pch 0 space
    typedef Lagrange<0, Scalar, Discontinuous> lag_0_basis_type;
    typedef FunctionSpace<mesh_type, bases<lag_0_basis_type>> lag_0_space_type;
    typedef std::shared_ptr<lag_0_space_type> lag_0_space_ptrtype;
    typedef typename lag_0_space_type::element_type lag_0_element_type;

    //! Pchv space
    typedef Lagrange<1, Vectorial> lag_v_basis_type;
    typedef FunctionSpace<mesh_type, bases<lag_v_basis_type>> lag_v_space_type;
    typedef std::shared_ptr<lag_v_space_type> lag_v_space_ptrtype;
    typedef typename lag_v_space_type::element_type lag_v_element_type;

    typedef FunctionSpace<mesh_type, bases<curl_basis_type,lag_basis_type>> comp_space_type;
    typedef std::shared_ptr<comp_space_type> comp_space_ptrtype;
    typedef typename comp_space_type::element_type comp_element_type;

    //! Preconditioners
    typedef PreconditionerBlockMS<comp_space_type> preconditioner_type;
    typedef std::shared_ptr<preconditioner_type> preconditioner_ptrtype;

    //! The exporter factory
    typedef Exporter<mesh_type> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    //! Backends factory
    typedef Backend<double> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef backend_type::solve_return_type solve_ret_type;

    public:

    /// Init the geometry with a circle/sphere from radius and characteristic length
    ///     \param radius   Circle or sphere radius.
    ///     \param h        Mesh size.
    TestPrecAFP( ) 
    {
        auto M_mesh = loadMesh(_mesh=new mesh_type);
        auto Xh = comp_space_type::New(M_mesh); // curl x lag
        auto Mh = lag_0_space_type::New(M_mesh); // lag_0
        auto MVh= lag_v_space_type::New(M_mesh); // lag_v
        auto M_mu_r = Mh->element();
        auto model = ModelProperties(Environment::expand(soption("myModel")));
        /*
         * we construct:
         * B = mu_r mu_0 H = mu H
         * we evaluate this relation for H=1
         * thus : B(mu_r,mu_0,H) = mu
         */
        for(auto it : model.materials() ){
            std::string key = "Materials.";
            key += marker(it);
            key += ".B";
            auto curve = expr(model.getEntry(key),{"H","mu_r"},{expr("1."),expr(soption("functions.m"))});
            curve.setParameterValues(model.parameters().toParameterValues());
            M_mu_r += vf::project(_space=Mh,
                                _range=markedelements(M_mesh,marker(it)),
                                _expr=curve/model.parameters().toParameterValues()["mu_0"]);
        }
        auto f_M_a = expr<DIM,1,6>(soption("functions.a"));
        auto c_M_a = expr<DIM,1,6>(soption("functions.c"));
        auto rhs = expr<DIM,1,6>(soption("functions.j"));
        std::pair<std::string, double> p;
        p.first = "m";
        p.second = std::stod(soption("functions.m"));
        rhs.setParameterValues(p);
        auto M_rhs = vf::project(_space=MVh, _range=elements(M_mesh), _expr=(rhs));

        auto U = Xh->element();
        auto V = Xh->element();
        auto u = U.template element<0>();
        auto v = V.template element<0>();
        auto phi = U.template element<1>();
        auto psi = V.template element<1>();

        auto f2 = form2(_test=Xh,_trial=Xh);
        auto f1 = form1(_test=Xh);

        preconditioner_ptrtype M_prec;
       
        map_vector_field<DIM,1,2> m_w_u {model.boundaryConditions().getVectorFields<DIM>("u","Weakdir")};
        map_scalar_field<2> m_w_phi {model.boundaryConditions().getScalarFields<2>("phi","Weakdir")};

        map_vector_field<DIM,1,2> m_dirichlet {model.boundaryConditions().getVectorFields<DIM>("u","Dirichlet")};
        map_scalar_field<2> m_dirichlet_phi {model.boundaryConditions().getScalarFields<2>("phi","Dirichlet")};
        
        f1 = integrate(_range= elements(M_mesh),
                       _expr = inner(idv(M_rhs),id(v)));    // rhs
        f2 = integrate(_range= elements(M_mesh),
                       _expr = 
                         inner(trans(id(v)),gradt(phi)) // grad(phi)
                       + inner(trans(idt(u)),grad(psi)) // div(u) = 0
                       + (1./idv(M_mu_r))*(trans(curlt_op(u))*curl_op(v)) // curl curl
                       );
        for(auto const & it : m_w_u)
        {
            LOG(INFO) << "Setting (weak) " << it.second << " on " << it.first << std::endl;
            f1 += integrate(_range=markedfaces(M_mesh,it.first), _expr=
                -(cst(1.)/idv(M_mu_r))*trans(curl_op(v))*cross(N(),it.second) 
                + doption("penaldir")/(idv(M_mu_r)*hFace())*inner(cross(it.second,N()),cross(id(v),N())) );
            f2 += integrate(_range=markedfaces(M_mesh,it.first), 
                _expr=-(1/idv(M_mu_r))*trans(curlt_op(u))*(cross(N(),id(v)) )
                - (1/idv(M_mu_r))*trans(curl_op(v))*(cross(N(),idt(u)) )
                + doption("penaldir")/(hFace()*idv(M_mu_r))*inner(cross(idt(u),N()),cross(id(v),N())) );
        }
        for(auto const & it : m_w_phi)
        {
                LOG(INFO) << "Setting (weak) " << it.second << " on " << it.first << std::endl;
                f2 += integrate(_range=markedfaces(M_mesh,it.first), 
                    _expr=doption("penaldir")/(hFace())*inner(idt(phi),id(psi)) );
        }
        for(auto const & it : m_dirichlet)
        {
            LOG(INFO) << "Setting " << it.second << " on " << it.first << std::endl;
            f2 += on(_range=markedfaces(M_mesh,it.first),
                     _rhs=f1,
                     _element=u,
                     _expr=it.second);
        }
        for(auto const & it : m_dirichlet_phi)
        {
            LOG(INFO) << "Setting " << it.second << " on " << it.first << std::endl;
            f2 += on(_range=markedfaces(M_mesh,it.first),
                     _rhs=f1,
                     _element=phi,
                     _expr=it.second);
        }

        solve_ret_type ret;
        if(soption("ms.pc-type") == "ms" ){
            // auto M_prec = blockms(
            //    _space = Xh,
            //    _space2 = Mh,
            //    _matrix = f2.matrixPtr(),
            //    _bc = model.boundaryConditions());

            M_prec = std::make_shared<PreconditionerBlockMS<comp_space_type>>(Xh, 
                                                                                model,
                                                                                "ms",
                                                                                f2.matrixPtr(), 0.1);

            //M_prec->update(f2.matrixPtr(),M_mu_r);
            tic();
            ret = f2.solveb(_rhs=f1,
                      _solution=U,
                      _backend=backend(_name="ms"),
                      _prec=M_prec);
            toc("Inverse",FLAGS_v>0);
        }else{
            tic();
            ret = f2.solveb(_rhs=f1,
                      _solution=U,
                      _backend=backend(_name="ms"));
            toc("Inverse",FLAGS_v>0);
        }
#if 1
        Environment::saveTimers(boption("saveTimers")); 
        auto e21 = normL2(_range=elements(M_mesh), _expr=(f_M_a-idv(u)));
        auto e22 = normL2(_range=elements(M_mesh), _expr=f_M_a);
        auto e21_curl = integrate(_range=elements(M_mesh),
                                  _expr = inner(idv(u)-f_M_a,
                                                idv(u)-f_M_a)
                                  + inner(curlv(u)-c_M_a,
                                          curlv(u)-c_M_a)
                                 ).evaluate()(0,0);
        auto e22_curl = integrate(_range=elements(M_mesh),
                                  _expr = inner(idv(u),idv(u))
                                  + inner(curlv(u),curlv(u))
                                 ).evaluate()(0,0);
        
        auto ecurl1 = normL2(_range=elements(M_mesh), _expr=(c_M_a-curlv(u)));
        auto ecurl2 = normL2(_range=elements(M_mesh), _expr=c_M_a);
        
        Feel::cout << "Error\t" << e21 << "\t" << e21/e22 << "\t" << e21_curl << "\t" << e21_curl/e22_curl << "\t" << ecurl1 << "\t" << ecurl1/ecurl2 << std::endl;
#else
        auto nnzVec = f2.matrixPtr()->graph()->nNz();
        int nnz = std::accumulate(nnzVec.begin(),nnzVec.end(),0);
        M_prec->iterMinMaxMean();
        if(Environment::worldComm().globalRank()==0)
            /*
             * to print * -> ./myCode --options >> sameResFile.txt
             * + grep #TestCase > parseIt
             * #TestCase
             * #ksp-pc_ksp11-pc11_ksp11.1-pc11.1_ksp11.2-pc11.2_ksp22-pc22
             * hSize
             * nProc
             * nDof
             * nDof(space1)
             * nDof(space2)
             * iter total
             * iter block11 : min max mean
             * iter block22 : min max mean
             * timer iter total
             * timer iter block11 : min max mean
             * timer iter block22 : min max mean
             * mu : min max mean
             */
        std::Feel::cout 
            << ioption("test-case") << "\t"
            << soption("ms.ksp-type") << "-" << soption("ms.pc-type") << "_"
            << soption("blockms.11.ksp-type") << "-" << soption("blockms.11.pc-type") << "_"
            << soption("blockms.11.1.ksp-type") << "-" << soption("blockms.11.1.pc-type") << "_"
            << soption("blockms.11.2.ksp-type") << "-" << soption("blockms.11.2.pc-type") << "_"
            << soption("blockms.22.ksp-type") << "-" << soption("blockms.22.pc-type") << "\t"
            << doption("gmsh.hsize") << "\t"
            //<< nProc << "\t"
            << Xh->nDof() << "\t"
            << Xh->template functionSpace<0>()->nDof() << "\t"
            << Xh->template functionSpace<1>()->nDof() << "\t"
            << ret.nIterations() << "\t"
            << M_prec->printMinMaxMean(0,0) << "\t"
            << M_prec->printMinMaxMean(0,1) << "\t"
            << M_prec->printMinMaxMean(0,2) << "\t"
            << M_prec->printMinMaxMean(1,0) << "\t"
            << M_prec->printMinMaxMean(1,1) << "\t"
            << M_prec->printMinMaxMean(1,2) << "\t"
            << e21 << "\t"
            << e21/e22 << "\n"; 
            //std::Feel::cout << "RES\t"
            //    << Xh->nDof() << "\t"
            //    << nnz << "\t"
            //    << soption("functions.m") << "\t"
            //    << e21 << "\t"
            //    << e21/e22
#if 0
            //    << "\t"
            //    << e21_curl << "\t"
            //    << e21_curl/e22_curl
#endif
            //    << std::endl;
        /* report */
        if ( Environment::worldComm().isMasterRank() && boption("generateMD") ){
            time_t now = std::time(0);
            tm *ltm = localtime(&now);
            std::ostringstream stringStream;
            stringStream << 1900+ltm->tm_year<<"_"<<ltm->tm_mon<<"_"<<ltm->tm_mday<<"-"<<ltm->tm_hour<<ltm->tm_min<<ltm->tm_sec<<"_"<<soption("title")<<".md";
            std::ofstream outputFile( stringStream.str() );
            if( outputFile )
            {
                outputFile
                        << "---\n"
                        << "title: \""<< soption("title") << "\"\n"
                        << "date: " << stringStream.str() << "\n"
                        << "categories: simu\n"
                        << "--- \n\n";
                outputFile << "#Physique" << std::endl;
                model.saveMD(outputFile);

                outputFile << "##Physique specifique" << std::endl;
                outputFile << "| Variable | value | " << std::endl;
                outputFile << "|---|---|" << std::endl;
                outputFile << "| mu | "    << soption("functions.m") << "| " << std::endl;
                outputFile << "| Rhs | "   << soption("functions.j") << "|" << std::endl;
                outputFile << "| Exact | " << soption("functions.a") << "|" << std::endl;

                outputFile << "#Numerics" << std::endl;

                outputFile << "##Mesh" << std::endl;
                M_mesh->saveMD(outputFile);

                outputFile << "##Spaces" << std::endl;
                outputFile << "|qDim|" << Xh->qDim()      << "|" << Xh->template functionSpace<1>()->qDim()      << "|" << Xh->template functionSpace<1>()->qDim()      << "|" << std::endl;
                outputFile << "|---|---|---|---|" << std::endl;
                outputFile << "|BasisName|"<< Xh->basisName() << "|" << Xh->template functionSpace<0>()->basisName() << "|" << Xh->template functionSpace<1>()->basisName() << "|" << std::endl;
                outputFile << "|nDof|" << Xh->nDof()      << "|"<< Xh->template functionSpace<0>()->nDof()      << "|"<< Xh->template functionSpace<1>()->nDof()      << "|" << std::endl;
                outputFile << "|nLocaldof|" << Xh->nLocalDof() << "|" << Xh->template functionSpace<0>()->nLocalDof() << "|" << Xh->template functionSpace<1>()->nLocalDof() << "|" << std::endl;
                outputFile << "|nPerComponent|" << Xh->nDofPerComponent() << "|" << Xh->template functionSpace<0>()->nDofPerComponent() << "|" << Xh->template functionSpace<1>()->nDofPerComponent() << "|" << std::endl;

                outputFile << "##Solvers" << std::endl;

                outputFile << "| x | ms | blocksms.11 | blockms.22 |" << std::endl;
                outputFile << "|---|---|---|---| " << std::endl;
                outputFile << "|**ksp-type** |  " << soption("ms.ksp-type") << "| " << soption("blockms.11.ksp-type") << "| " << soption("blockms.22.ksp-type") << "|" << std::endl;
                outputFile << "|**pc-type**  |  " << soption("ms.pc-type")  << "| " << soption("blockms.11.pc-type")  << "| " << soption("blockms.22.pc-type")  << "|" << std::endl;

                if(soption("ms.pc-type") == "blockms" ){
                    outputFile << "|**Matrix**  |  " << nnz << "| "; M_prec->printMatSize(1,outputFile); outputFile << "| "; M_prec->printMatSize(2,outputFile);outputFile  << "|" << std::endl;
                    outputFile << "|**nb Iter**  |  " << ret.nIterations() << "| "; M_prec->printIter(1,outputFile); outputFile << "| "; M_prec->printIter(2,outputFile);outputFile  << "|" << std::endl;
                }else{
                    outputFile << "|**Matrix**  |  " << nnz << "| 0 | 0 |" << std::endl;
                    outputFile << "|**nb Iter**  |  " << ret.nIterations() << "| 0 | 0 |" << std::endl;
                }

                outputFile << "##Timers" << std::endl;
                Environment::saveTimersMD(outputFile);
            }
            else
            {
               std::cerr << "Failure opening " << stringStream.str() << '\n';
            }
        }
        /* end of report */
#endif
        // export
        if(boption("exporter.export")){
            auto exact = vf::project(_space=Xh->template functionSpace<0>(),_range=elements(M_mesh), _expr=f_M_a);
            auto diff  = vf::project(_space=Xh->template functionSpace<0>(),_range=elements(M_mesh), _expr=idv(u)-f_M_a);
            auto ex = exporter(_mesh=M_mesh);
            ex->add("relativPermeability",M_mu_r  );
            ex->add("rhs"                ,M_rhs  );
            ex->add("potential"          ,u  );
            ex->add("potential_e"        ,exact  );
            ex->add("potential_d"        ,diff  );
            ex->add("lagrange_multiplier",phi);
            ex->save();
        }
    }

private:
    mesh_ptrtype M_mesh;
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( precAFP )

BOOST_AUTO_TEST_CASE( test )
{
    TestPrecAFP<FEELPP_DIM> test;
}

//// Test 3D
//BOOST_AUTO_TEST_CASE( test_3d )
//{
//    TestPrecAFP<3> test;
//}

BOOST_AUTO_TEST_SUITE_END()

#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _about=makeAbout(),
                           _desc=makeOptions() );

    TestPrecAFP<FEELPP_DIM> t_afp;

    return 0;
}
#endif
