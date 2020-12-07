#ifndef FEELPP_MODELS_STABILIZATIONGLSPARAMETERBASE_HPP
#define FEELPP_MODELS_STABILIZATIONGLSPARAMETERBASE_HPP 1

#include <feel/feeldiscr/pdh.hpp>

namespace Feel {

namespace FeelModels
{


template<typename MeshType>
class StabilizationGLSParameterBase
{
public :
    using mesh_t = MeshType;
    using mesh_ptr_t = std::shared_ptr<mesh_t>;

    using space_tau_t = Pdh_type<mesh_t,0>;
    using space_tau_ptr_t = std::shared_ptr<space_tau_t>;
    using element_tau_t = typename space_tau_t::element_type;
    using element_tau_ptr_t = typename space_tau_t::element_ptrtype;

    StabilizationGLSParameterBase( mesh_ptr_t const& mesh, std::string const& prefix, po::variables_map const& vm )
        :
        M_mesh( mesh ),
        M_method( soption(_name="method",_prefix=prefix,_vm=vm ) ),
        M_hSizeMethod( soption(_name="hsize.method",_prefix=prefix,_vm=vm ) ),
        M_penalLambdaK( doption( _name="eigenvalue.penal-lambdaK",_prefix=prefix,_vm=vm ) )
        {
            CHECK( M_method == "eigenvalue" || M_method == "doubly-asymptotic-approximation" )  << "invalid method";
        }

    virtual ~StabilizationGLSParameterBase() = default;
    
    virtual void init() = 0;

    double hSize( size_type eltId ) const
        {
            auto itFindElt = M_hSizeValues.find( eltId );
            CHECK(itFindElt != M_hSizeValues.end() ) << "no elt id " << eltId << " found in hSize";
            return itFindElt->second;
        }
    double sqrtLambdaK(size_type eltId ) const
        {
            auto itFindElt = M_sqrtLambdaKValues.find( eltId );
            CHECK(itFindElt != M_sqrtLambdaKValues.end() ) << "no elt id " << eltId << " found in sqrtLambdaK";
            return itFindElt->second;
        }
    double lambdaK( size_type eltId ) const
        {
            auto itFindElt = M_lambdaKValues.find( eltId );
            CHECK(itFindElt != M_lambdaKValues.end() ) << "no elt id " << eltId << " found in lambdaK";
            return itFindElt->second;
        }
    double mK( size_type eltId ) const
        {
            auto itFindElt = M_mKValues.find( eltId );
            CHECK(itFindElt != M_mKValues.end() ) << "no elt id " << eltId << " found in mK";
            return itFindElt->second;
        }

    mesh_ptr_t const& mesh() const { return M_mesh; }
    std::string const& method() const { return M_method; }
    std::string const& hSizeMethod() const { return M_hSizeMethod; }
    double penalLambdaK() const { return M_penalLambdaK; }
    void setPenalLambdaK( double val ) { M_penalLambdaK = val; }

    element_tau_ptr_t const& fieldTauPtr() const { return M_fieldTau; }
    element_tau_t const& fieldTau() const { return *M_fieldTau; }

protected :
    mesh_ptr_t M_mesh;
    std::unordered_map<size_type,double> M_hSizeValues, M_lambdaKValues, M_sqrtLambdaKValues, M_mKValues;
    std::string M_method;
    std::string M_hSizeMethod;
    double M_penalLambdaK;

    space_tau_ptr_t M_spaceTau;
    element_tau_ptr_t M_fieldTau;
    //element_tau_ptr_t M_fieldDelta;
};

}  // namespace FeelModels

} // namespace Feel

#endif // FEELPP_MODELS_STABILIZATIONGLSPARAMETERBASE_HPP
