/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_HEAT3D_HPP
#define FEELPP_HEAT3D_HPP 1

#include <feel/options.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeHeat3dOptions();
FEELPP_EXPORT AboutData
makeHeat3dAbout( std::string const& str = "heat3d" );

struct FEELPP_EXPORT Heat3dConfig
{
    typedef Mesh<Simplex<3>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;
};

class FEELPP_EXPORT Heat3d : public ModelCrbBase<ParameterSpace<>, Heat3dConfig::space_type >
{
    typedef ModelCrbBase<ParameterSpace<>, Heat3dConfig::space_type > super_type;

public:
    Heat3d();

    void initBetaQ();

    super_type::betaq_type computeBetaQ( parameter_type const& mu );


    void updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const;

    //! initialisation of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);


    // parameter_type crbParametersFromUserParameters( feelpp4fastsim::UserParameters const& userParam ) const;
    // void updateFieldsExported( SourceCrbInterface * pvsource, element_type & feField, vtkSmartPointer<vtkUnstructuredGrid> vtkOutput );


};

}

#endif /* FEELPP_HEAT3D_HPP */
