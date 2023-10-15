/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

#ifndef FEELPP_OPUSHEAT_HPP
#define FEELPP_OPUSHEAT_HPP 1

#include <boost/timer.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelts/bdf.hpp>


namespace Feel
{

FEELPP_EXPORT po::options_description
makeOpusHeatOptions();
FEELPP_EXPORT AboutData
makeOpusHeatAbout( std::string const& str = "opusheat" );

template<bool IsStationary>
struct FEELPP_EXPORT OpusHeatConfig
{
    typedef Mesh<Simplex<2>> mesh_type;
    typedef Pch_type<mesh_type,1> space_type;
    static const int options = (IsStationary)? TimeIndependent : TimeDependent;
};


template<bool IsStationary>
class FEELPP_EXPORT OpusHeat : public ModelCrbBase< ParameterSpace<>, typename OpusHeatConfig<IsStationary>::space_type, OpusHeatConfig<IsStationary>::options >
{
    typedef ModelCrbBase<ParameterSpace<>, typename OpusHeatConfig<IsStationary>::space_type , OpusHeatConfig<IsStationary>::options > super_type;
public:

    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;

    using element_type = typename super_type::element_type;
    using element_ptrtype = std::shared_ptr<element_type>;
    using parameter_type = typename super_type::parameter_type;
    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
//    using rbspace_type = typename super_type::rbspace_type;
//    using rbspace_ptrtype = typename super_type::rbspace_ptrtype;

    using super_type::computeBetaQm;

   // using bdf_type = std::variant< Bdf<space_type>, Bdf<rbspace_type> >;

    using bdf_type = Bdf<space_type>;
    using bdf_ptrtype = std::shared_ptr<bdf_type>;

    OpusHeat();

    void initModel();

    typename super_type::betaq_type
    computeBetaQ( parameter_type const& mu , double time , bool only_terms_time_dependent=false );
    typename super_type::betaq_type
    computeBetaQ( parameter_type const& mu );

    double output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false );


    bdf_ptrtype bdfModel(){ return M_bdf; }

    void initializationField( element_ptrtype& initial_field,parameter_type const& mu );

    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir );
    void updateSpecificityModel( boost::property_tree::ptree & ptree ) const;


  private :
    void initDataStructureAffineDecomposition();
    void assembleData();

    void updateBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent );

    template<typename ModelType>
    typename super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaAq, this->M_betaFq );
        }

    template<typename ModelType>
    typename super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu , double time , bool only_terms_time_dependent=false,
                       typename std::enable_if< ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,time,only_terms_time_dependent);
            return boost::make_tuple( this->M_betaMq, this->M_betaAq, this->M_betaFq );
        }


    template<typename ModelType>
    typename super_type::betaq_type
    computeBetaQ_impl( parameter_type const& mu, typename std::enable_if< !ModelType::is_time_dependent >::type* = nullptr )
        {
            this->updateBetaQ_impl( mu,0,false);
            return boost::make_tuple( this->M_betaAq, this->M_betaFq );
        }

    template<typename ModelType>
    typename super_type::betaq_type
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

    std::map<std::string,double> M_measureMarkedSurface;
    bdf_ptrtype M_bdf;

};

}

#endif /* FEELPP_OPUSHEAT_HPP */
