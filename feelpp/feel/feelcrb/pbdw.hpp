#ifndef _FEELPP_PBDW_HPP
#define _FEELPP_PBDW_HPP 1
#include <Eigen/Core>
namespace Feel
{

using Eigen::MatrixXd;
template<typename CRB_Model, typename SensorMap >
class PBDW
{
    typedf CRB_Model::element_type element_type;

public:

    PBDW( int U_max_dim, bool do_offline = true):
        M_M_max(U_max_dim),
        M_do_offline(do_offline),
        M_M(1),
        {}
    ~PBDW()
        {}

    static std::vector<double> get_measurement( double time, std::vector<std::string> capteurs_name, int numberOfMeasure);

    void offline()
    {
        // construction of Z_N
        if ( offlineStep()==false )
        {
            auto conv = CRB_Model.offline();
        }

        auto N = CRB_model.dimension();
        auto Z_N  = CRB_Model.wn();
        auto res = CRB_Model.orthonormalize(N, Z_N, N);
        // construction of U_M;
        std::array<element_type, M_M_max> UM;
        auto a = CRB_Model.algebraicInnerProduct();
        i = 0;
        for( auto const& sensor : SensorMap )
        {
            element_type q;
            lm = sensor;
            a.solve(_rhs=lm, _solution=q);
            UM(i) = q;
            i ++;
        }

        MatrixXd A(M_M,M_M);

        for(int i=0; i<M_M; i++)
        {
          for(int j=0; j<M_M; j++)

              A(i,j) = CRB_Model.scalarProduct(UM[j], UM[i]);
        }

        MatrixXd B(M_M,N);

        for(int i=0; i<M_M; i++)
        {
          for(int j=0; j<N; j++)

            B(i,j) = CRB_Model.scalarProduct(ZN[j], UM[i]);
        }


        MatrixXd K(M_M+N,M_M+N);

        K.setZero();
        K.block(0,0,M_M,M_M) = A;
        K.block(0,M_M,M_M,N) = B;
        K.block(M_M,0,N,M_M) = B.transpose();

        // Save K
        MatrixXd I_U(M_M,1);
        for(int i=0; i<M_M; i++)
          I_U(i) = CRB_Model.out(UM[i]);

        for(int i=0; i<N; i++)
          I_Z(i) = CRB_Model.out(ZN[i]);

        //Save I_U and I_Z
    }

    element_type online(std::array<double,M_M> y_obs)
    {

        MatrixXd Y(M_M+N,1);
        Y.setZero();
        for(int i=0; i<M_M; i++)
          Y(i,0)=y_obs[i];

        MatrixXd S(M_M+N,1);
        S = K.inverse()*Y; // load matrix K

        auto eta = S.block(0,0,M_M,0);
        auto z = S.block(M_M,0,N,0);

        auto ouput = eta * I_U + z * I_Z; // load I_U and I_Z

    }

protected

    int M_M_max, M_M;
    bool M_do_offline;
};



}






#endif /* _FEELPP_PBDW_HPP */
