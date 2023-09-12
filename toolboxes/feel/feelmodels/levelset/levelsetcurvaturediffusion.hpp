/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <metivet@math.unistra.fr>
 Date: 2018-11-05

 Copyright (C) 2018 Université de Strasbourg

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
 \file levelsetcurvaturediffusion.hpp
 \author Thibaut Metivet <metivet@math.unistra.fr>
 \date 2018-11-05
 */
#ifndef _LEVELSETCURVATUREDIFFUSION_HPP
#define _LEVELSETCURVATUREDIFFUSION_HPP 1

namespace Feel {
namespace FeelModels {

template< typename FunctionSpaceType >
class LevelSetCurvatureDiffusion
{
    typedef LevelSetCurvatureDiffusion< FunctionSpaceType > self_type;

public:
    typedef FunctionSpaceType functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    static inline const uint16_type Dim = functionspace_type::nDim;
    static inline const uint16_type Order = functionspace_type::basis_type::nOrder;

    typedef typename functionspace_type::element_type element_type;
    typedef typename functionspace_type::element_ptrtype element_ptrtype;

    typedef typename functionspace_type::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef elements_reference_wrapper_t<typename functionspace_type::mesh_type> range_elements_type;
    typedef faces_reference_wrapper_t<typename functionspace_type::mesh_type> range_faces_type;

    enum class TimeDiscretisation { CRANK_NICOLSON, EULER1 };

public:
    LevelSetCurvatureDiffusion( functionspace_ptrtype space, std::string const& prefix );
    virtual ~LevelSetCurvatureDiffusion() {}

    std::string const& prefix() const { return M_prefix; }
 
    functionspace_ptrtype const& functionSpace() const { return M_functionSpace; }
    range_elements_type rangeMeshElements() const { return M_rangeMeshElements; }
    range_faces_type rangeMeshBoundaryFaces() const { return M_rangeMeshBoundaryFaces; }

    double timeStep() const { return M_timeStep; }
    void setTimeStep( double dt );

    TimeDiscretisation timeDiscretisation() const { return M_timeDiscretisation; }
    void setTimeDiscretisation( TimeDiscretisation d );

    element_type solveGalpha( element_type const& phi ) const;
    element_type solveGalpha( element_ptrtype const& phi ) const { return this->solveGalpha( *phi ); }
    element_type solveGbeta( element_type const& phi ) const;
    element_type solveGbeta( element_ptrtype const& phi ) const { return this->solveGbeta( *phi ); }

    element_type curvatureOrder1( element_type const& phi ) const;
    element_type curvatureOrder1( element_ptrtype const& phi ) const { return this->curvatureOrder1( *phi ); }
    element_type curvatureOrder2( element_type const& phi ) const;
    element_type curvatureOrder2( element_ptrtype const& phi ) const { return this->curvatureOrder2( *phi ); }

    element_type willmore( element_type const& phi ) const { return this->willmoreImpl( phi, mpl::int_<Dim>() ); }
    element_type willmore( element_ptrtype const& phi ) const { return this->willmore( *phi ); }

private:
    void initCurvatureDiffusionAlpha();
    void initCurvatureDiffusionBeta();

    element_type willmoreImpl( element_type const& phi, mpl::int_<1> /*Dim*/ ) const {}
    element_type willmoreImpl( element_type const& phi, mpl::int_<2> /*Dim*/ ) const;
    element_type willmoreImpl( element_type const& phi, mpl::int_<3> /*Dim*/ ) const;

private:
    std::string M_prefix;

    functionspace_ptrtype M_functionSpace;
    range_elements_type M_rangeMeshElements;
    range_faces_type M_rangeMeshBoundaryFaces;

    // alpha = sqrt(2)
    static const double M_alpha;
    // beta = 1/sqrt(2)
    static const double M_beta;

    double M_timeStep;

    TimeDiscretisation M_timeDiscretisation;

    // Note: the backends are duplicated to avoid a weird crash with PETSC/MUMPS
    // which seems to hold information regarding the matrix, and crashes when solving
    // both the alpha and beta systems with the same backend. 
    // This should be investigated further...
    backend_ptrtype M_backendAlpha;
    backend_ptrtype M_backendBeta;
    matrix_ptrtype M_curvatureDiffusion_alphaDt;
    matrix_ptrtype M_curvatureDiffusion_betaDt;
    mutable vector_ptrtype M_vector;

};

template<typename FunctionSpaceType>
const double LevelSetCurvatureDiffusion<FunctionSpaceType>::M_alpha = std::sqrt(2.);
template<typename FunctionSpaceType>
const double LevelSetCurvatureDiffusion<FunctionSpaceType>::M_beta = 1./std::sqrt(2.);

template<typename FunctionSpaceType>
LevelSetCurvatureDiffusion<FunctionSpaceType>::LevelSetCurvatureDiffusion( functionspace_ptrtype space, std::string const& prefix )
    :
        M_prefix( prefix ),
        M_functionSpace( space ),
        M_rangeMeshElements( elements( space->mesh() ) ),
        M_rangeMeshBoundaryFaces( boundaryfaces( space->mesh() ) ),
        M_backendAlpha( backend_type::build( soption( _name="backend", _prefix=prefix ), prefix, space->worldCommPtr() ) ),
        M_backendBeta( backend_type::build( soption( _name="backend", _prefix=prefix ), prefix, space->worldCommPtr() ) )
{
    if( Environment::vm( _name="time-step", _prefix=this->prefix() ).defaulted() )
        M_timeStep = std::pow( this->functionSpace()->mesh()->hAverage()/Order, 2);
    else
        M_timeStep = doption( _name="time-step", _prefix=this->prefix() );

    std::string timeDiscretisation = soption( _name="time-discretisation", _prefix=this->prefix() );
    if( timeDiscretisation == "crank-nicolson" )
        M_timeDiscretisation = TimeDiscretisation::CRANK_NICOLSON;
    else if( timeDiscretisation == "euler1" )
        M_timeDiscretisation = TimeDiscretisation::EULER1;
    else
        CHECK( false ) << timeDiscretisation << " is not in the list of available time discretisations\n";

    this->initCurvatureDiffusionAlpha();
    this->initCurvatureDiffusionBeta();
    M_vector = M_backendAlpha->newVector( this->functionSpace() );
}

template<typename FunctionSpaceType>
void
LevelSetCurvatureDiffusion<FunctionSpaceType>::setTimeStep( double dt )
{
    M_timeStep = dt;
    this->initCurvatureDiffusionAlpha();
    this->initCurvatureDiffusionBeta();
}

template<typename FunctionSpaceType>
void
LevelSetCurvatureDiffusion<FunctionSpaceType>::setTimeDiscretisation( TimeDiscretisation d )
{
    M_timeDiscretisation = d;
    this->initCurvatureDiffusionAlpha();
    this->initCurvatureDiffusionBeta();
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::solveGalpha( element_type const& phi ) const
{
    auto linearForm = form1( _test=this->functionSpace(), _vector=M_vector );
    /* Galpha */
    switch( M_timeDiscretisation )
    {
        case TimeDiscretisation::CRANK_NICOLSON:
        {
            // Crank-Nicolson scheme
            linearForm = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idv(phi) * id(phi) / (M_alpha*M_timeStep) - 0.5 * inner(gradv(phi), grad(phi))
                    );
        } break;
        case TimeDiscretisation::EULER1:
        {
            // Euler1 scheme
            linearForm = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idv(phi) * id(phi) / (M_alpha*M_timeStep) 
                    );
        } break;
    }
    // Explicit Neumann BC
    linearForm += integrate(
            _range=this->rangeMeshBoundaryFaces(),
            _expr=id(phi) * gradv(phi) * vf::N()
            );
    // Solve
    auto Galpha = this->functionSpace()->element();
    M_backendAlpha->solve( _matrix=M_curvatureDiffusion_alphaDt, _rhs=M_vector, _solution=Galpha );

    return Galpha;
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::solveGbeta( element_type const& phi ) const
{
    auto linearForm = form1( _test=this->functionSpace(), _vector=M_vector );
    /* Gbeta */
    switch( M_timeDiscretisation )
    {
        case TimeDiscretisation::CRANK_NICOLSON:
        {
            // Crank-Nicolson scheme
            linearForm = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idv(phi) * id(phi) / (M_beta*M_timeStep) - 0.5 * inner(gradv(phi), grad(phi))
                    );
        } break;
        case TimeDiscretisation::EULER1:
        {
            // Euler1 scheme
            linearForm = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idv(phi) * id(phi) / (M_beta*M_timeStep)
                    );
        } break;
    }
    // Explicit Neumann BC
    linearForm += integrate(
            _range=this->rangeMeshBoundaryFaces(),
            _expr=id(phi) * gradv(phi) * vf::N()
            );
    // Solve
    auto Gbeta = this->functionSpace()->element();
    M_backendBeta->solve( _matrix=M_curvatureDiffusion_betaDt, _rhs=M_vector, _solution=Gbeta );

    return Gbeta;
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::curvatureOrder1( element_type const& phi ) const
{
    // Solve Gbeta
    auto Gbeta = this->solveGbeta( phi );

    return vf::project( _space=this->functionSpace(), _range=this->rangeMeshElements(), 
            _expr=(idv(Gbeta)-idv(phi))/(M_beta*M_timeStep) 
            );
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::curvatureOrder2( element_type const& phi ) const
{
    /* Galpha */
    auto Galpha = this->solveGalpha( phi );
    /* Gbeta */
    auto Gbeta = this->solveGbeta( phi );

    return vf::project( _space=this->functionSpace(), _range=this->rangeMeshElements(), 
            _expr=(-idv(Galpha) + 4*idv(Gbeta) - 3*idv(phi))/(M_alpha*M_timeStep) 
            );
}

template<typename FunctionSpaceType>
void
LevelSetCurvatureDiffusion<FunctionSpaceType>::initCurvatureDiffusionAlpha()
{
    if( !M_curvatureDiffusion_alphaDt )
        M_curvatureDiffusion_alphaDt = M_backendAlpha->newMatrix( _trial=this->functionSpace(), _test=this->functionSpace() );

    auto bilinearForm_alphaDt = form2( 
            _trial=this->functionSpace(),
            _test=this->functionSpace(),
            _matrix=M_curvatureDiffusion_alphaDt
            );

    auto phi = this->functionSpace()->element();
    switch( M_timeDiscretisation )
    {
        case TimeDiscretisation::CRANK_NICOLSON:
        {
            // Heat equation with Crank-Nicolson scheme
            bilinearForm_alphaDt = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idt(phi) * id(phi) / (M_alpha*M_timeStep) + 0.5 * inner(gradt(phi), grad(phi))
                    );
        } break;
        case TimeDiscretisation::EULER1:
        {
            // Heat equation with Euler1 scheme
            bilinearForm_alphaDt = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idt(phi) * id(phi) / (M_alpha*M_timeStep) + inner(gradt(phi), grad(phi))
                    );
        } break;
    }
}

template<typename FunctionSpaceType>
void
LevelSetCurvatureDiffusion<FunctionSpaceType>::initCurvatureDiffusionBeta()
{
    if( !M_curvatureDiffusion_betaDt )
        M_curvatureDiffusion_betaDt = M_backendBeta->newMatrix( _trial=this->functionSpace(), _test=this->functionSpace() );

    auto bilinearForm_betaDt = form2( 
            _trial=this->functionSpace(),
            _test=this->functionSpace(),
            _matrix=M_curvatureDiffusion_betaDt
            );

    auto phi = this->functionSpace()->element();
    switch( M_timeDiscretisation )
    {
        case TimeDiscretisation::CRANK_NICOLSON:
        {
            // Heat equation with Crank-Nicolson scheme
            bilinearForm_betaDt = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idt(phi) * id(phi) / (M_beta*M_timeStep) + 0.5 * inner(gradt(phi), grad(phi))
                    );
        } break;
        case TimeDiscretisation::EULER1:
        {
            // Heat equation with Euler1 scheme
            bilinearForm_betaDt = integrate(
                    _range=this->rangeMeshElements(),
                    _expr=idt(phi) * id(phi) / (M_beta*M_timeStep) + inner(gradt(phi), grad(phi))
                    );
        } break;
    }
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::willmoreImpl( element_type const& phi, mpl::int_<2> /*Dim*/ ) const
{
    /* Galpha */
    auto Galpha = this->solveGalpha( phi );
    /* Gbeta */
    auto Gbeta = this->solveGbeta( phi );

    return vf::project( _space=this->functionSpace(), _range=this->rangeMeshElements(), 
            _expr=(idv(phi) + idv(Galpha) - 2*idv(Gbeta))/(M_timeStep*M_timeStep) 
            );
}

template<typename FunctionSpaceType>
typename LevelSetCurvatureDiffusion<FunctionSpaceType>::element_type
LevelSetCurvatureDiffusion<FunctionSpaceType>::willmoreImpl( element_type const& phi, mpl::int_<3> /*Dim*/ ) const
{
    //TODO
    return this->functionSpace()->element();
}

     
} // namespace FeelModels
} // namespace Feel

#endif

