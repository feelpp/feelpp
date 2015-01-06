/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Cecile Daversin  <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-12-07

  Copyright (C) 2011 UJF
  Copyright (C) 2011 CNRS

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file grepl_fem.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2014-12-22
 */
#include <feel/feel.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description greploptions( "test grepl_fem" );
    greploptions.add_options()
        ( "mu1", po::value<double>()->default_value( 0.01 ), "param 1" )
        ( "mu2", po::value<double>()->default_value( 0.01 ), "param 2" )
        ( "gamma", po::value<double>()->default_value( 10 ), "weak dirichlet" )
        ( "tol", po::value<double>()->default_value( 1e-12 ), "picard tolerance" )
        ( "weakdir", po::value<bool>()->default_value( true ), "picard tolerance" )
        ;
    return greploptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "grepl_fem" ,
                     "grepl_fem" ,
                     "0.1",
                     "solve problem of grepl's benchmark (for rb use)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

class GreplFem
    :
public Application
{
    typedef Application super;

public:

    typedef GreplFem self_type;

    //! numerical type is double
    typedef double value_type;
    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<2,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    //typedef bases<Lagrange<1,Scalar> > basis_type; //Lagrange scalar space
    typedef bases<Lagrange<2,Scalar> > basis_type; //Lagrange scalar space
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /**
     * Constructor
     */
    GreplFem()
        :
        super(),
        //M_backend( backend_type::build( this->vm() ) ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        this->changeRepository( boost::format( "/grepl_fem/%1%-%2%" )
                                % doption(_name="mu1")
                                % doption(_name="mu2")
                                );
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _desc=domain( _name = "grepl_fem_mesh",
                                             _shape = "hypercube",
                                             //_usenames = true,
                                             _dim = 2,
                                             _h=2e-2,
                                             _xmin=0,_xmax=1,
                                             _ymin=0,_ymax=1 ) );
    }

    void run();
    void updateJacobian( vector_ptrtype const& X, sparse_matrix_ptrtype & J);
    void updateResidual(vector_ptrtype const& X,vector_ptrtype & R);
    element_type solveNewton();
    element_type solvePicard();

    private:
    //! linear algebra backend
    //backend_ptrtype M_backend;
    //! exporter factory
    export_ptrtype exporter;
    //! mesh
    mesh_ptrtype mesh;

}; //GreplFem

void
GreplFem::run()
{
    export_ptrtype exporter( export_type::New( this->vm(), "grepl_fem" ) );
    exporter->step( 0 )->setMesh( mesh );
    std::cout << "newton..." << std::endl;
    //call solveNewton and export
    auto solNewton = solveNewton();
    std::cout << "min newton = " << solNewton.min() << std::endl;
    std::cout << "max newton = " << solNewton.max() << std::endl;
    std::cout << "output = " << integrate( elements(mesh), idv(solNewton) ).evaluate()(0,0) << std::endl;

    exporter->step( 0 )->add( "newton", solNewton );
    std::cout << "picard..." << std::endl;
    //call solvePicard and export
    auto solPicard = solvePicard();
    std::cout << "min picard = " << solPicard.min() << std::endl;
    std::cout << "max picard = " << solPicard.max() << std::endl;
    std::cout << "output = " << integrate( elements(mesh), idv(solPicard) ).evaluate()(0,0) << std::endl;

    exporter->step( 0 )->add( "picard", solPicard );
    std::cout << "end!" << std::endl;

    exporter->save();
}

void
GreplFem::updateJacobian( vector_ptrtype const& X, sparse_matrix_ptrtype & J)
{
    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    double gamma = option(_name="gamma").template as<double>();
    double mu1=doption(_name="mu1");
    double mu2=doption(_name="mu2");

    auto g = exp( mu2*idv(u) );
    //auto proj_g = vf::project(_space=Xh, _expr=g);

    J->zero();
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range= elements( mesh ),_expr = gradt(u)*trans(grad(v)) );
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = boundaryfaces(mesh),
                                                          _expr = gamma*idt(u)*id(v)/hFace()
                                                          - (gradt(u)*vf::N())*id(v)
                                                          - (grad(v)*vf::N())*idt(u) );

    //form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = elements(mesh), _expr = mu1*idv(proj_g)*idt(u)*id(v) );
    form2( _test=Xh, _trial=Xh, _matrix=J ) += integrate( _range = elements(mesh), _expr = mu1*g*idt(u)*id(v) );

    J->close();
}

void
GreplFem::updateResidual(vector_ptrtype const& X,vector_ptrtype & R)
{
    auto Xh = space_type::New( mesh );
    double gamma = option(_name="gamma").template as<double>();
    double mu1=doption(_name="mu1");
    double mu2=doption(_name="mu2");

    auto u = Xh->element();
    u=*X;
    auto v = Xh->element(); //test

    auto g = exp( mu2*idv(u) );

    R->zero();
    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ), _expr = gradv(u)*trans(grad(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range = boundaryfaces( mesh ),
                   _expr = gamma*idv(u)*id(v)/hFace()
                   - (gradv(u)*vf::N())*id(v) +
                   - (grad(v)*vf::N())*idv(u) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   _expr= mu1/mu2*( g * id(v) ) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   _expr= -mu1/mu2*( id(v)) );

    form1( _test=Xh, _vector=R ) +=
        integrate( _range= elements( mesh ),
                   //_expr=-100*sin(2*M_PI*Px())*cos(2*M_PI*Py()) * id(v) );
                   _expr=-100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );

    R->close();

}

GreplFem::element_type
GreplFem::solveNewton()
{
    auto Xh = space_type::New( mesh );
    auto J = backend()->newMatrix( Xh, Xh);
    auto R = backend()->newVector( Xh );
    backend()->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian,
                                                   boost::ref( *this ), _1, _2);
    backend()->nlSolver()->residual = boost::bind( &self_type::updateResidual,
                                                   boost::ref( *this ), _1, _2);

    auto solution = Xh->element();
    backend()->nlSolve(_jacobian=J, _solution=solution, _residual=R);

    return solution;

}

GreplFem::element_type
GreplFem::solvePicard()
{
    auto Xh = space_type::New( mesh );
    double gamma = option(_name="gamma").template as<double>();
    double mu1=doption(_name="mu1");
    double mu2=doption(_name="mu2");
    double tol = doption(_name="tol");
    bool weakdir = boption(_name="weakdir");
    int nb_iter_max = 100;

    auto solution = Xh->element();
    double error=0.0;
    int iter=0;
    do{
        auto u=Xh->element();
        auto v=Xh->element();
        auto solution_old = solution; //store previous iteration

        auto exprg = exp( mu2*idv(solution) );

        auto a = form2( _test=Xh, _trial=Xh );
        a = integrate( _range= elements( mesh ), _expr = gradt(u)*trans(grad(v)) );
        a += integrate( _range= elements( mesh ), _expr= mu1/mu2*( exprg*id(v) - id(v) ));

        if( weakdir )
        {
            a += integrate( _range = boundaryfaces( mesh ),
                            _expr = gamma*idt(u)*id(v)/hFace()
                            - (gradt(u)*vf::N())*id(v)
                            - (grad(v)*vf::N())*idt(u) );
        }

        auto rhs = form1( _test=Xh );
        rhs = integrate( _range= elements( mesh ), _expr = 100*sin(2*M_PI*Px())*sin(2*M_PI*Py()) * id(v) );

        if( !weakdir )
            a += on(_range=boundaryfaces(mesh), _rhs=rhs, _element=solution, _expr=cst(0.) );

        a.solve( _rhs=rhs, _solution=solution );

        error = normL2( _range=elements(mesh), _expr=idv(solution_old) - idv(solution));
        std::cout << "iter = " << iter << ", error = " << error << std::endl;
        iter++;

    }while(error > tol && iter < nb_iter_max);

    return solution;
}

int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv, makeOptions() );
    GreplFem app;
    app.run();
}
