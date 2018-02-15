/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_LINEARELASTICITY3D_HPP
#define FEELPP_LINEARELASTICITY3D_HPP 1

#include <feel/options.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeLinearElasticity3dOptions();
FEELPP_EXPORT AboutData
makeLinearElasticity3dAbout( std::string const& str = "linearelasticity3d"/*"linear-elasticity-3d"*/ );

struct FEELPP_EXPORT LinearElasticity3dConfig
{
    typedef Mesh<Simplex<3>> mesh_type;
    typedef Pchv_type<mesh_type,1> space_type;
};

class FEELPP_EXPORT LinearElasticity3d : public ModelCrbBase<ParameterSpace<>, LinearElasticity3dConfig::space_type >
{
    typedef ModelCrbBase<ParameterSpace<>, LinearElasticity3dConfig::space_type > super_type;
public:

    LinearElasticity3d();

    void initBetaQ();

    super_type::betaq_type computeBetaQ( parameter_type const& mu );

    //! initialisation of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

};


}

#endif /* FEELPP_LINEARELASTICITY3D_HPP */
