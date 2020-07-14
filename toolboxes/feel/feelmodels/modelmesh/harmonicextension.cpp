/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Copyright (C) 2014 Université Joseph Fourier (Grenoble I)

  This library Is free software; you can redistribute it and/or
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
   \file harmonicextension.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2014-05-16
*/

#include <feel/feelmodels/modelmesh/harmonicextension.hpp>

#include <feel/feeldiscr/stencil.hpp>

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/unary.hpp>
#include <feel/feelvf/one.hpp>
#include <feel/feelvf/cst.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/trans.hpp>
#include <feel/feelvf/inner.hpp>
#include <feel/feelvf/trace.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/val.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/on.hpp>

#include <feel/feelmodels/modelcore/log.hpp>


namespace Feel
{
namespace FeelModels
{

template< typename MeshType, int Order >
HarmonicExtension<MeshType,Order>::HarmonicExtension( mesh_ptrtype mesh, backend_ptrtype const& backend, std::string prefix,
                                                      worldcomm_ptr_t const& worldcomm,
                                                      bool useGhostEltFromExtendedStencil,
                                                      ModelBaseRepository const& modelRep )
    :
    super_type( prefix, worldcomm,"",modelRep ),
    ModelBase( prefix, worldcomm,"",modelRep ),
    M_backend( backend ),
    M_mesh( mesh ),
    M_useAdaptPenal( option(_prefix=this->prefix(),_name="use_adaptive_penalisation").template as<bool>() )
{
    M_Xh = space_type::New(_mesh=this->mesh(),_worldscomm=makeWorldsComm(1,worldcomm),
                           _extended_doftable=std::vector<bool>(1,useGhostEltFromExtendedStencil) );
    M_displacement = M_Xh->elementPtr();
    M_dispImposedOnBoundary = M_Xh->elementPtr();
    M_XhP0 = space_P0_type::New( _mesh=this->mesh(),_worldscomm=makeWorldsComm(1,worldcomm) );
}

template< typename MeshType, int Order >
HarmonicExtension<MeshType,Order>::HarmonicExtension( space_ptrtype space, backend_ptrtype const& backend, std::string prefix,
                                                      ModelBaseRepository const& modelRep )
    :
    super_type( prefix, space->worldCommPtr(),"",modelRep ),
    ModelBase( prefix, space->worldCommPtr(),"",modelRep ),
    M_backend( backend ),
    M_mesh( space->mesh() ),
    M_Xh( space ),
    M_useAdaptPenal( option(_prefix=this->prefix(),_name="use_adaptive_penalisation").template as<bool>() )
{
    M_displacement = M_Xh->elementPtr();
    M_dispImposedOnBoundary = M_Xh->elementPtr();
    M_XhP0 = space_P0_type::New( _mesh=this->mesh(),_worldscomm=makeWorldsComm(1,this->worldCommPtr()) );
}


template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::backend_ptrtype const&
HarmonicExtension<MeshType,Order>::backend() const { return M_backend; }

template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::mesh_ptrtype const&
HarmonicExtension<MeshType,Order>::mesh() const { return M_mesh; }

template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::space_ptrtype const&
HarmonicExtension<MeshType,Order>::functionSpace() const { return M_Xh; }

template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::element_ptrtype const&
HarmonicExtension<MeshType,Order>::displacement() const { return M_displacement; }

template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::element_ptrtype const&
HarmonicExtension<MeshType,Order>::dispImposedOnBoundary() const { return M_dispImposedOnBoundary; }

template< typename MeshType, int Order >
typename HarmonicExtension<MeshType,Order>::flagSet_type const&
HarmonicExtension<MeshType,Order>::flagSet() const { return M_flagSet; }

template< typename MeshType, int Order >
std::vector<flag_type> const&
HarmonicExtension<MeshType,Order>::flagSet(std::string key) const
{
    CHECK( M_flagSet.find(key) != M_flagSet.end() ) << "the flag type " << key << " is unknown \n";
    return M_flagSet.find(key)->second;
}
template< typename MeshType, int Order >
flag_type
HarmonicExtension<MeshType,Order>::flagSet(std::string key, int k) const
{
    CHECK( M_flagSet.find(key) != M_flagSet.end() ) << "the flag type " << key << " is unknown \n";
    CHECK( M_flagSet.find(key)->second.size() > k ) << "the key " << k << " must be <  " <<  M_flagSet.find(key)->second.size() << "\n";
    return M_flagSet.find(key)->second.at(k);
}

template< typename MeshType, int Order >
void
HarmonicExtension<MeshType,Order>::setflagSet( flagSet_type const & fl )
{
    M_flagSet=fl;
}

#if 0 // code for implement near null space
namespace detail
{
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<2> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
#if 1
    auto mode3 = space->element( vec(Py(),-Px()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
#else
    NullSpace<double> userNullSpace( { mode1,mode2 } );
#endif
    return userNullSpace;
}
template <typename SpaceType>
NullSpace<double> getNullSpace( SpaceType const& space, mpl::int_<3> /**/ )
{
    auto mode1 = space->element( oneX() );
    auto mode2 = space->element( oneY() );
    auto mode3 = space->element( oneZ() );
#if 1
    auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
    auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
    auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
    NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
#else
    NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
#endif
    return userNullSpace;
}
} // namespace detail
#endif

template< typename MeshType, int Order >
void
HarmonicExtension<MeshType,Order>::init()
{
    this->log(this->prefix(),"init", "start");
    // mesh not move so not rebuild cst part
    this->setRebuildCstPartInLinearSystem(false);
#if 1
    auto graph = stencil( _test=M_Xh, _trial=M_Xh )->graph();
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend(),
                                                         graph, graph->mapRow().indexSplit() ) );
#else
    // A bug to fix with this code ( probably block pattern )
    M_algebraicFactory.reset( new model_algebraic_factory_type( this->shared_from_this(),this->backend() ) );
    M_algebraicFactory->initFromFunctionSpace( this->functionSpace() );
#endif

#if 0 // code for implement near null space
    if ( true )
    {
        NullSpace<double> userNullSpace = detail::getNullSpace(this->functionSpace(), mpl::int_<mesh_type::nDim>() ) ;
        //if ( boption(_name="use-null-space",_prefix=this->prefix() ) )
        //M_algebraicFactory->attachNullSpace( userNullSpace );
        M_algebraicFactory->attachNearNullSpace( userNullSpace );
    }
#endif

    M_vectorSolution = this->backend()->newVector( this->functionSpace() );

    this->log(this->prefix(),"init", "finish");
}

template< typename MeshType, int Order >
std::shared_ptr<std::ostringstream>
HarmonicExtension<MeshType,Order>::getInfo() const
{
    std::shared_ptr<std::ostringstream> _ostr( new std::ostringstream() );

    *_ostr << "\n   Harmonic Extension "
           << "\n     -- type : " << ((M_useAdaptPenal)?"yes":"no")
        ;

    return _ostr;
}


template< typename MeshType, int Order >
void
HarmonicExtension<MeshType,Order>::updateLinearPDE( DataUpdateLinear & data ) const
{
    //const vector_ptrtype& X = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();

    auto u = this->displacement();
    auto v = this->displacement();

    if ( buildCstPart )
    {
        if ( M_useAdaptPenal )
        {
            //std::cout << " harmonic : use_adaptive_penalisation \n";
            /* Calculate penalization term due to Masud and Hughes \tau_e = \frac{1-V_{min}/V_{max}}{V_e/V_{max}} */
            auto XhP0 = M_XhP0;
            auto Volume = integrate( elements(this->mesh()), meas() ).broken( XhP0 );
            double Vmin = Volume.min();
            double Vmax = Volume.max();
            auto tau = ((1.0 - Vmin/Vmax)/( idv(Volume)/Vmax ));

            form2( _test=M_Xh, _trial=M_Xh, _matrix=A ) +=
                integrate( _range=elements(this->mesh()),
                           _expr=val(1+tau)*trace( trans(gradt(u))*grad(v) ) );
#if 0
            for ( uint16_type i=0; i < M_flagSet["free"].size(); ++i )
                form1( _test=M_Xh, _vector=F ) +=
                    integrate(_range=markedfaces(this->mesh(), this->flagSet("free",i)),
                              _expr=val(1+tau)*trans(eye<Dim,Dim>()/*gradv(P())*/*N())*id(v) );
#endif
            //form2( _test=fspaceLow, _trial=fspaceLow, _matrix=harmonicLow ) +=
            //    integrate( _range=boundaryfaces(reference_mesh), _expr= -val(1+tau)*trans((gradt(u)*N()))*id(v) );
        }
        else
        {
            //std::cout << " standart harmonic\n";
            form2( _test=M_Xh, _trial=M_Xh, _matrix=A ) +=
                integrate( _range=elements(this->mesh()),
                           _expr=trace( trans(gradt(u))*grad(v) ) );
#if 0
            for ( uint16_type i=0; i < M_flagSet["free"].size(); ++i )
                form1( _test=M_Xh, _vector=F ) +=
                    integrate(_range=markedfaces(this->mesh(), this->flagSet("free",i)),
                              _expr=trans(eye<Dim,Dim>()/*gradt(u)*/*N())*id(v) );
#endif
            //form2( _test=fspaceLow, _trial=fspaceLow, _matrix=harmonicLow ) +=
            //    integrate( _range=boundaryfaces(reference_mesh), _expr= -trans((gradt(u)*N()))*id(v) );
        }
    } // if ( buildCstPart )

    if ( !buildCstPart )
    {
        for ( uint16_type i=0; i < this->flagSet("moving").size(); ++i )
        {
            //if (this->worldComm().isMasterRank() ) std::cout << "M_flagSet[moving][i]" << M_flagSet["moving"][i]<< std::endl;
            form2( _test=M_Xh, _trial=M_Xh, _matrix=A ) +=
                on( _range=markedfaces( this->mesh(), this->flagSet("moving",i) ),
                    _element=*u,
                    _rhs=F,
                    _expr=idv( this->dispImposedOnBoundary() ) );
        }

        for ( uint16_type i=0; i < this->flagSet("fixed").size(); ++i )
        {
            //if (this->worldComm().isMasterRank() ) std::cout << "M_flagSet[fixed][i]" << M_flagSet["fixed"][i]<< std::endl;
            form2( _test=M_Xh, _trial=M_Xh, _matrix=A ) +=
                on( _range=markedfaces(this->mesh(), this->flagSet("fixed",i) ),
                    _element=*u,
                    _rhs=F,
                    _expr=0*one() );
        }
    }

}

template< typename MeshType, int Order >
void
HarmonicExtension<MeshType,Order>::solve()
{
    this->log(this->prefix(),"solve", "start");
    // assemble and solve linear system
    M_algebraicFactory->solveLinear(M_vectorSolution);
    // copy algebraic vector into an finite element approximation
    *M_displacement = *M_vectorSolution;

    this->log(this->prefix(),"solve", "finish");
}


template class HarmonicExtension< Mesh<Simplex<2,1> >, 1 >;
template class HarmonicExtension< Mesh<Simplex<3,1> >, 1 >;

} // namespace FeelModels
} // namespace Feel

