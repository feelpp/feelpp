#ifndef FEELPP_STOKES_HPP
#define FEELPP_STOKES_HPP 1

#include <feel/options.hpp>
#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feelmor/crbsaddlepoint.hpp>
#include <feel/feelmor/crbmodelsaddlepoint.hpp>
#include <feel/feeldiscr/thch.hpp>


namespace Feel
{

FEELPP_EXPORT AboutData
    makeStokesDeimAbout( std::string const& str = "stokesdeim" );

FEELPP_EXPORT po::options_description
makeStokesDeimOptions()
{
    po::options_description stokesDeimOptions( "StokesDeim options" );
    return stokesDeimOptions;
}

struct FEELPP_EXPORT StokesDeimConfig
{
    typedef Mesh<Simplex<2>> mesh_type;
    typedef THch_type<1,mesh_type> space_type;
};

class FEELPP_EXPORT StokesDeim :
    public ModelCrbBase<ParameterSpace<>, StokesDeimConfig::space_type, UseBlock >
{
    typedef StokesDeim self_type;
    typedef ModelCrbBase<ParameterSpace<>, StokesDeimConfig::space_type, UseBlock > super_type;

 public :
    typedef typename StokesDeimConfig::space_type space_type;
    typedef boost::tuple<beta_vector_type,  std::vector<beta_vector_type> > beta_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::displacement_space_type displacement_space_type;
    typedef typename super_type::displacement_space_ptrtype displacement_space_ptrtype;
    typedef typename super_type::displacement_field_ptrtype displacement_field_ptrtype;
    using super_type::computeBetaQm;

    StokesDeim();
    void initModel() override;
    beta_type computeBetaQm( parameter_type const& mu ) override;
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false) override;
    sparse_matrix_ptrtype assembleForMDEIM( parameter_type const& mu, int const& tag ) override;
    displacement_field_ptrtype meshDisplacementField( parameter_type const& mu ) override;

  private:
    displacement_space_ptrtype M_Dh;
    displacement_field_ptrtype mapping;
};


} // namespace Feel

#endif
