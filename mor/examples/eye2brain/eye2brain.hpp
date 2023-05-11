/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_EYE2BRAIN_HPP
#define FEELPP_EYE2BRAIN_HPP 1

#define ORDER 1

#include <feel/options.hpp>

#include <feel/feelmor/modelcrbbase.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeEye2BrainOptions();
FEELPP_EXPORT AboutData
makeEye2BrainAbout( std::string const& str = "eye2brain" );

template<int Order, int Dim>
struct FEELPP_EXPORT Eye2BrainConfig
{
    typedef double value_type;
    typedef Mesh<Simplex<Dim>> mesh_type;
    typedef bases< Lagrange<Order, Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    // typedef Pch_type<mesh_type, 2> space_type;
};

template<int Dim>
class FEELPP_EXPORT Eye2Brain : public ModelCrbBase<ParameterSpace<>, typename Eye2BrainConfig<ORDER, Dim>::space_type >
{
    typedef ModelCrbBase<ParameterSpace<>, typename Eye2BrainConfig<ORDER, Dim>::space_type > super_type;
public:
    typedef CRBResults::parameter_type parameter_type;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::value_type value_type;

    Eye2Brain();

    void initBetaQ();

    typename super_type::betaq_type computeBetaQ( parameter_type const& mu );

    void updateSpecificityModelPropertyTree( boost::property_tree::ptree & ptree ) const;

    //! initialisation of the model
    void initModel();

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    double output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);


    // parameter_type crbParametersFromUserParameters( feelpp4fastsim::UserParameters const& userParam ) const;
    // void updateFieldsExported( SourceCrbInterface * pvsource, element_type & feField, vtkSmartPointer<vtkUnstructuredGrid> vtkOutput );


};

}

#endif /* FEELPP_EYE2BRAIN_HPP */
