/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_THERMALBUILDING_HPP
#define FEELPP_THERMALBUILDING_HPP 1

#include <feel/options.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeThermalBuildingOptions();
FEELPP_EXPORT AboutData
makeThermalBuildingAbout( std::string const& str = "thermalbuilding" );

struct FEELPP_EXPORT ThermalBuildingConfig
{
    typedef Mesh<Simplex<3>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;
    static const int options = TimeIndependent;

};

class FEELPP_EXPORT ThermalBuilding : public ModelCrbBase<ParameterSpace<>, ThermalBuildingConfig::space_type >
{
    typedef ModelCrbBase<ParameterSpace<>, ThermalBuildingConfig::space_type, ThermalBuildingConfig::options > super_type;
public:

    ThermalBuilding();

    //! initialisation of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);

    void initBetaQ();

    super_type::betaq_type
    computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent=false );
    super_type::betaq_type
    computeBetaQ( parameter_type const& mu );


private :
    void initDataStructureAffineDecomposition();
    void assembleData();

private:
    int M_Qa;
    std::vector< int > M_Ql;
    int M_Nl;
};

}

#endif /* FEELPP_THERMALBUILDING_HPP */
