/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Thibaut Metivet <metivet@math.unistra.fr>
 Date: 2019-03-04

 Copyright (C) 2019 Université de Strasbourg

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
 \file levelsetparticleinjector.hpp
 \author Thibaut Metivet <metivet@math.unistra.fr>
 \date 2019-03-04
 */
#ifndef _LEVELSETPARTICLEINJECTOR_HPP
#define _LEVELSETPARTICLEINJECTOR_HPP 1

#include<feel/feelmodels/levelset/levelsetparticleshapes.hpp>

namespace Feel {
namespace FeelModels {

template< typename LevelsetType >
class LevelSetParticleInjector
{
public:
    typedef LevelSetParticleInjector<LevelsetType> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef LevelsetType levelset_type;
    typedef std::shared_ptr<levelset_type> levelset_ptrtype;

    typedef typename levelset_type::value_type value_type;
    //--------------------------------------------------------------------//
    // Mesh
    static inline const uint16_type nDim = levelset_type::nDim;
    typedef typename levelset_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    // Range types
    typedef typename MeshTraits<mesh_type>::element_reference_wrapper_const_iterator element_reference_wrapper_const_iterator;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_type elements_reference_wrapper_type;
    typedef typename MeshTraits<mesh_type>::elements_reference_wrapper_ptrtype elements_reference_wrapper_ptrtype;
    typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;
    // Levelset element
    typedef typename levelset_type::space_levelset_type functionspace_type;
    typedef typename levelset_type::space_levelset_ptrtype functionspace_ptrtype;
    typedef typename functionspace_type::element_type element_levelset_type;
    typedef typename functionspace_type::element_ptrtype element_levelset_ptrtype;
    // Levelset particle shapes
    typedef LevelSetParticleShapes<functionspace_type> levelsetparticleshapes_type;
    typedef std::shared_ptr<levelsetparticleshapes_type> levelsetparticleshapes_ptrtype; 
    // Particle injection method
    enum class ParticleInjectionMethod {
        FIXED_POSITION
    };
    static const std::map<std::string, ParticleInjectionMethod> ParticleInjectionMethodIdMap;
    // Particles
    struct Particle {
        LevelSetShapeType shape;
        parameter_map parameters;
    };

    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//
    //--------------------------------------------------------------------//

    //--------------------------------------------------------------------//
    // Constructor
    LevelSetParticleInjector( levelset_ptrtype const& ls, range_elements_type const& range );

    //--------------------------------------------------------------------//
    // Function spaces
    functionspace_ptrtype const& functionSpaceLevelset() const { return M_levelset->functionSpace(); }
    functionspace_ptrtype const& functionSpaceInjector() const { return M_spaceInjector; }

    //--------------------------------------------------------------------//
    // Parameters
    range_elements_type const& rangeInjectorElements() const { return M_rangeInjectorElements; }
    //void setInjectorMarkers( boost::any const& markers ) { M_injectorMarkers = markers; }

    LevelSetShapeType particleShape( std::string const& id ) const { return M_particles.at(id).shape; }
    void setParticleShape( std::string const& id, LevelSetShapeType shape ) { M_particles[id].shape = shape; }
    void setParticleShape( std::string const& id, std::string const& shape );

    parameter_map const& particleParams( std::string const& id ) const { return M_particles.at(id).parameters; }
    void setParticleParams( std::string const& id, parameter_map const& params ) { M_particles[id].parameters = params; }

    void addParticle( std::string const& id, LevelSetShapeType shape, parameter_map const& params );
    void addParticle( std::string const& id, std::string const& shape, parameter_map const& params );

    ParticleInjectionMethod particleInjectionMethod() const { return M_particleInjectionMethod; }
    void setParticleInjectionMethod( ParticleInjectionMethod method ) { M_particleInjectionMethod = method; }
    void setParticleInjectionMethod( std::string const& method );

    uint16_type nParticles() const { return M_nParticles; }
    void setNParticles( uint16_type n ) { M_nParticles = n; }

    double particleDensityThreshold() const { return M_particleDensityThreshold; }
    void setParticleDensityThreshold( double t ) { M_particleDensityThreshold = t; }

    parameter_map const& fixedParticleParams() const { return M_fixedParticleParams; }
    void setFixedParticleParams( parameter_map const& params ) { M_fixedParticleParams = params; }

    //--------------------------------------------------------------------//
    // Particle shapes
    levelsetparticleshapes_ptrtype const& levelsetParticleShapes() const { return M_levelsetParticleShapes; }
    void setLevelsetParticleShapes( levelsetparticleshapes_ptrtype const& s ) { M_levelsetParticleShapes = s; }

    //--------------------------------------------------------------------//
    // Inject
    element_levelset_type inject( element_levelset_type const& ls ) const;

private:
    //--------------------------------------------------------------------//
    // Levelset
    levelset_ptrtype M_levelset;
    //--------------------------------------------------------------------//
    // Injector function space
    range_elements_type M_rangeInjectorElements;
    functionspace_ptrtype M_spaceInjector;
    //--------------------------------------------------------------------//
    // Particle shapes creator
    levelsetparticleshapes_ptrtype M_levelsetParticleShapes;
    //--------------------------------------------------------------------//
    // Parameters
    std::map<std::string, Particle> M_particles;
    ParticleInjectionMethod M_particleInjectionMethod;
    uint16_type M_nParticles;
    double M_particleDensityThreshold;
    parameter_map M_fixedParticleParams;

};

#define LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS \
    template< typename LevelsetType > \
    /**/
#define LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE \
    LevelSetParticleInjector< LevelsetType > \
    /**/

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
const std::map<std::string, typename LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::ParticleInjectionMethod>
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::ParticleInjectionMethodIdMap = {
    { "fixed_position", ParticleInjectionMethod::FIXED_POSITION }
};

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::LevelSetParticleInjector( levelset_ptrtype const& ls, range_elements_type const& range )
    : M_levelset( ls ),
      M_rangeInjectorElements( range ),
      M_particles(),
      M_particleInjectionMethod( ParticleInjectionMethod::FIXED_POSITION ),
      M_nParticles( 0 ),
      M_particleDensityThreshold( 1e-6 )
{
    //M_spaceInjector = functionspace_type::New( 
            //_mesh=this->functionSpaceLevelset()->mesh(), 
            //_worldscomm=this->functionSpaceLevelset()->worldsComm(), 
            //_range=markedelements(this->functionSpaceLevelset()->mesh(), M_injectorMarkers)
            //);
    // Temporary hack
    M_spaceInjector = this->functionSpaceLevelset();
    M_levelsetParticleShapes = std::make_shared<levelsetparticleshapes_type>( M_spaceInjector );
}

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::setParticleShape( std::string const& id, std::string const& shape )
{
    CHECK( LevelSetShapeTypeIdMap.count( shape ) ) << shape << " is not in the list of supported particle shapes\n";
    this->setParticleShape( id, LevelSetShapeTypeIdMap.at( shape ) );
}

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::addParticle( std::string const& id, LevelSetShapeType shape, parameter_map const& params )
{
    M_particles[id] = { shape, params };
}

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::addParticle( std::string const& id, std::string const& shape, parameter_map const& params )
{
    CHECK( LevelSetShapeTypeIdMap.count( shape ) ) << shape << " is not in the list of supported particle shapes\n";
    this->addParticle( id, LevelSetShapeTypeIdMap.at(shape), params );
}

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
void
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::setParticleInjectionMethod( std::string const& method )
{
    CHECK( ParticleInjectionMethodIdMap.count( method ) ) << method << " is not in the list of supported particle injection methods\n";
    this->setParticleInjectionMethod( ParticleInjectionMethodIdMap.at( method ) );
}

LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
typename LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::element_levelset_type
LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE::inject( element_levelset_type const& phi ) const
{
    auto const& mesh = phi.mesh();

    double injectionZoneVolume = integrate(
            _range=this->rangeInjectorElements(),
            _expr=cst(1.)
            ).evaluate()(0,0);
    double particleDensity = integrate(
            _range=this->rangeInjectorElements(),
            _expr=(1-idv(M_levelset->heaviside()))
            ).evaluate()(0,0) / injectionZoneVolume;

    auto phiInject = phi;

    if( particleDensity <= this->particleDensityThreshold() )
    {
        element_levelset_type phiParticles = vf::project(
                _space=this->functionSpaceInjector(),
                _range=this->functionSpaceInjector()->template rangeElements<0>(),
                _expr=idv(phi)
                );
        switch( M_particleInjectionMethod )
        {
            case ParticleInjectionMethod::FIXED_POSITION:
            {
                for( auto const& particle: M_particles )
                {
                    auto const newPart = this->levelsetParticleShapes()->create( particle.second.shape, particle.second.parameters, true );
                    phiParticles = vf::project(
                            _space=this->functionSpaceInjector(),
                            _range=this->functionSpaceInjector()->template rangeElements<0>(),
                            _expr=vf::min( idv(phiParticles), idv(newPart) )
                            );
                }
            }
            break;
        }

        //phiInject.on( _range=rangeInjectorElements, _expr=idv(phiParticles) );
        // Temporary hack
        //phiInject = vf::project(
                //_space=phiInject.functionSpace(),
                //_range=phiInject.functionSpace()->template rangeElements<0>(),
                //_expr=vf::min( idv(phiParticles), idv(phiInject) )
                //);
        phiInject = phiParticles;
    }

    //return vf::project(
            //_space=this->functionSpaceLevelset(),
            //_range=elements(mesh),
            //_expr=vf::min( idv(phiParticles), idv(phi) )
            //);

    return phiInject;
}

#undef LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_DECLARATIONS
#undef LEVELSETPARTICLEINJECTOR_CLASS_TEMPLATE_TYPE

} // namespace FeelModels
} // namespace Feel

#endif // _LEVELSETPARTICLEINJECTOR_HPP

