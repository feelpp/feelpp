#ifndef _FEELPP_PBDW_HPP
#define _FEELPP_PBDW_HPP 1


namespace Feel
{
template<typename MeshType, typename FunctionSpaceType, typename ParameterSpaceType>
class ModelBase{};


template<typename ModelType, typename SensorMap, typename ScalarProduct, typename LinearForm >
class PBDW

{
    typedef ModelBase<SpaceType,typename ModelType::functionspace_type, typename ModelType::parameterspace_type> super;

public:

    typedef ModelType::functionspace_type::element_type element_type;
    typedef typename ModelType::functionspace_type functionspace_type;
    typedef typename ModelType::parameterspace_type  parameterspace_type;

    PBDW( int Z_max_dim , int U_max_dim, bool do_offline = true):
        M_M_max(U_max_dim),
        M_N_max(Z_max_dim),
        M_do_offline(do_offline),
        M_M(1),
        M_N(1)
        {}
    ~PBDW()
        {}


    void offline();

    element_type online( std::array<double,M_M> y_obs);

protected

    int M_M_max, M_N_max, M_M,M_N;
    bool M_do_offline;
};



}






#endif /* _FEELPP_PBDW_HPP */
