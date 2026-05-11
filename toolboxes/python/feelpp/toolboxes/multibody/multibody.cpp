//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

#include <feel/feelpython/pybind11/eigen.h>
#include <feel/feelpython/pybind11/pybind11.h>
#include <feel/feelpython/pybind11/stl.h>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/multibody/multibody.hpp>

namespace py = pybind11;
using namespace Feel;

template<typename ElementType>
py::object
castElementPtr( ElementType const& element )
{
    return py::cast( std::addressof( element ), py::return_value_policy::reference );
}

template<int nDim, int OrderGeo>
void
defMultibody( py::module &m )
{
    using namespace Feel;
    using namespace Feel::FeelModels;

    using toolbox_t = FeelModels::Multibody<Simplex<nDim, OrderGeo>>;
    using body_t = typename toolbox_t::body_type;

    std::string const body_class_name = fmt::format( "Body_{}DG{}", nDim, OrderGeo );
    py::class_<body_t>( m, body_class_name.c_str() )
        .def( "mass", &body_t::mass, "get the body mass" )
        .def( "massCenter", &body_t::massCenter, py::return_value_policy::reference_internal, "get the body mass center" )
        .def( "rigidTranslation", &body_t::rigidTranslation, py::return_value_policy::reference_internal, "get the body rigid translation" )
        .def( "rigidRotationAngles", &body_t::rigidRotationAngles, py::return_value_policy::reference_internal, "get the body rigid rotation angles" )
        .def( "momentOfInertia_bodyFrame", &body_t::momentOfInertia_bodyFrame, py::return_value_policy::reference_internal, "get the body-frame moment of inertia" )
        .def( "momentOfInertia_inertialFrame", &body_t::momentOfInertia_inertialFrame, "get the inertial-frame moment of inertia" )
        .def( "hasElasticDisplacement", &body_t::hasElasticDisplacement, "return true when elastic displacement is available" )
        .def( "hasElasticVelocity", &body_t::hasElasticVelocity, "return true when elastic velocity is available" )
        .def( "fieldDisplacement",
              []( body_t const& body ) { return castElementPtr( body.fieldDisplacement() ); },
              "get the body displacement field" )
        .def( "fieldElasticDisplacement",
              []( body_t const& body ) -> py::object
              {
                  if ( !body.hasElasticDisplacement() )
                      return py::none();
                  return castElementPtr( body.fieldElasticDisplacement() );
              },
              "get the body elastic displacement field, or None" )
        .def( "fieldElasticVelocity",
              []( body_t const& body ) -> py::object
              {
                  if ( !body.hasElasticVelocity() )
                      return py::none();
                  return castElementPtr( body.fieldElasticVelocity() );
              },
              "get the body elastic velocity field, or None" );

    std::string const pyclass_name = fmt::format( "Multibody_{}DG{}", nDim, OrderGeo );
    py::class_<toolbox_t, std::shared_ptr<toolbox_t>, ModelNumerical>( m, pyclass_name.c_str() )
        .def( py::init( []( std::string const& prefix, std::string const& keyword, py::object worldComm, ModelBaseRepository const& modelRep )
                        {
                            worldcomm_ptr_t wc = worldComm.is_none() ? Environment::worldCommPtr() : py::cast<worldcomm_ptr_t>( worldComm );
                            return new toolbox_t( prefix, keyword, wc, modelRep );
                        } ),
              py::arg( "prefix" ),
              py::arg( "keyword" ) = std::string( "multibody" ),
              py::arg( "worldComm" ) = py::none(),
              py::arg( "modelRep" ) = ModelBaseRepository(),
              "Initialize the multibody toolbox" )
        .def( "init", &toolbox_t::init, "initialize the multibody toolbox" )
        .def( "mesh", &toolbox_t::mesh, "get the mesh" )
        .def( "setMesh", &toolbox_t::setMesh, "set the mesh", py::arg( "mesh" ) )
        .def( "bodyNames",
              []( toolbox_t const& toolbox )
              {
                  std::vector<std::string> names;
                  names.reserve( toolbox.bodies().size() );
                  for ( auto const& [name, body] : toolbox.bodies() )
                      names.push_back( name );
                  return names;
              },
              "get the registered body names" )
        .def( "hasBody", &toolbox_t::hasBody, "return true when a body name is registered", py::arg( "name" ) )
        .def( "body",
              []( toolbox_t& toolbox, std::string const& name ) -> body_t&
              {
                  return toolbox.body( name );
              },
              py::return_value_policy::reference_internal,
              "get a body by name",
              py::arg( "name" ) );
}

PYBIND11_MODULE( _multibody, m )
{
    defMultibody<2, 1>( m );
    defMultibody<2, 2>( m );
    defMultibody<3, 1>( m );
    defMultibody<3, 2>( m );
}
