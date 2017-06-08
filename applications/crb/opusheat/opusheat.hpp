/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_OPUSHEAT_HPP
#define FEELPP_OPUSHEAT_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/feelts/bdf.hpp>


namespace Feel
{

po::options_description
makeOpusHeatOptions();
AboutData
makeOpusHeatAbout( std::string const& str = "opusheat" );

struct OpusHeatConfig
{
    typedef Mesh<Simplex<2>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;
    static const int options = TimeDependent;//TimeIndependent;//TimeDependent;//TimeDependent;//TimeIndependent;

    struct ParameterDefinition
    {
        static const uint16_type ParameterSpaceDimension = 6;
        typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    };

};


class OpusHeat : public ModelCrbBase< OpusHeatConfig::ParameterDefinition, OpusHeatConfig::space_type, OpusHeatConfig::options >
{
    typedef ModelCrbBase<OpusHeatConfig::ParameterDefinition, OpusHeatConfig::space_type , OpusHeatConfig::options > super_type;
public:

    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;

    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::parameter_type parameter_type;
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;

    using super_type::computeBetaQm;

    OpusHeat();

    void initModel();

    super_type::betaq_type
    computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent=false );
    super_type::betaq_type
    computeBetaQ( parameter_type const& mu );

    double output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false, bool export_outputs=false );


    bdf_ptrtype bdfModel(){ return M_bdf; }

    void initializationField( element_ptrtype& initial_field,parameter_type const& mu );
private :
    void initDataStructureAffineDecomposition();
    void assembleData();

    void updateBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent );

    template<typename ModelType>
    super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaAq, this->M_betaFq );
        }

    template<typename ModelType>
    super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaMq, this->M_betaAq, this->M_betaFq );
        }


    template<typename ModelType>
    super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu, typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,0,false);
            return boost::make_tuple( this->M_betaAq, this->M_betaFq );
        }

    template<typename ModelType>
    super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu, typename std::enable_if< ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,0,false);
            return boost::make_tuple( this->M_betaMq, this->M_betaAq, this->M_betaFq );
        }

private:

    int M_Qa;
    int M_Qm;
    std::vector< int > M_Ql;
    int M_Nl;

    double M_surface;
    bdf_ptrtype M_bdf;

};

}

#endif /* FEELPP_OPUSHEAT_HPP */
