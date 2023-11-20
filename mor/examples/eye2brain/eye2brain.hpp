/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_EYE2BRAIN_HPP
#define FEELPP_EYE2BRAIN_HPP 1

#define ORDER 2

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
    using value_type = double;
    using mesh_type = Mesh<Simplex<Dim>>;
    using basis_type = bases< Lagrange<Order, Scalar> >;
    using space_type = FunctionSpace<mesh_type, basis_type, value_type>;
    // typedef Pch_type<mesh_type, 2> space_type;
};

template<int Order, int Dim>
class FEELPP_EXPORT Eye2Brain : public ModelCrbBase<ParameterSpace<>, typename Eye2BrainConfig<Order, Dim>::space_type >
{
    using super_type = ModelCrbBase<ParameterSpace<>, typename Eye2BrainConfig<Order, Dim>::space_type >;
public:
    using parameter_type = typename super_type::parameter_type;
    using element_type = typename super_type::element_type;
    using value_type = typename super_type::value_type;
    using space_type = typename super_type::space_type;

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

private:
    const std::vector<std::vector<double>> m_coordinates = {
        //   O                  A                B                    B1                        C                  D               D1             F               G
            {-0.013597, 0, 0}, {-0.0126, 0, 0}, {-0.00979266, 0, 0}, {-0.009524,0.00194062,0}, {-0.00615819,0,0}, {0.0106, 0, 0}, {0.011, 0, 0}, {0.0115, 0, 0}, {0.0125, 0, 0}
        };
    const std::vector<std::string> m_outputNames = {"O", "A", "B", "B1", "C", "D", "D1", "F", "G"};
    double m_radius;
};

}

#endif /* FEELPP_EYE2BRAIN_HPP */
