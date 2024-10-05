#ifndef CRBMODELBLOCK_H
#define CRBMODELBLOCK_H

#include <feel/feelmor/crbmodel.hpp>

namespace Feel
{

template<typename ModelType>
class CRBModelBlock :
        public CRBModel<ModelType>
{
    typedef CRBModel<ModelType> super;
public:
    typedef ModelType model_type;
    typedef std::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;
    typedef typename model_type::parameter_type parameter_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::element_type element_type;

    template <int T>
    using subspace_type = typename space_type::template sub_functionspace<T>::type;
    template <int T>
    using subspace_ptrtype = std::shared_ptr<subspace_type<T>>;
    template<int T>
    using subelement_type = typename subspace_type<T>::element_type;

    static const int n_block = space_type::nSpaces;
    typedef typename mpl::range_c< int, 0, n_block > rangespace_type;

    CRBModelBlock( std::string modelName, crb::stage stage, int level = 0 ) :
        super ( modelName, stage, level ),
        M_block_initialized( false )
        {
            this->initBlocks();
        }

    CRBModelBlock( std::string modelName, model_ptrtype const& model , crb::stage stage, int level = 0 ) :
        super ( modelName, model, stage, level ),
        M_block_initialized( false )
        {
            this->initBlocks();
        }

    void initBlocks();

    using super::l2solve;
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f, int n_space )
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            M_backend_l2_vec[n_space]->solve( _matrix=M_inner_product_matrix_vec[n_space],
                                              _solution=u, _rhs=f );
        }

    void initBlockMatrix();
    void clearBlockMatrix();

    std::vector< std::vector< sparse_matrix_ptrtype >>
    AqmBlock( uint16_type n_space1, uint16_type n_space2 )
        {
            initBlockMatrix();
            return this->M_Aqm_block[n_space1][n_space2];
        }
    std::vector< std::vector< vector_ptrtype >>
    FqmBlock( uint16_type l, uint16_type n_space )
        {
            initBlockMatrix();
            return this->M_Fqm_block[l][n_space];
        }

    //!returns the scalar product of the vector x and vector y
    using super::scalarProduct;
    double scalarProduct( vector_type const& X, vector_type const& Y, int n_space )
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            return M_inner_product_matrix_vec[n_space]->energy( X, Y );
        }
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y, int n_space )
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            return M_inner_product_matrix_vec[n_space]->energy( X, Y );
        }

    virtual bool addSupremizerInSpace( int const& n_space ) const
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            return false;
        }

    virtual element_type supremizer( parameter_type const& mu, element_type const& U, int n_space )
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            return this->functionSpace()->element();
        }


protected:
    std::vector<sparse_matrix_ptrtype> M_inner_product_matrix_vec;
    std::vector< backend_ptrtype > M_backend_l2_vec;

    bool M_block_initialized;
    std::vector< std::vector< std::vector< std::vector<sparse_matrix_ptrtype>>>> M_Aqm_block;
    std::vector< std::vector< std::vector< std::vector<vector_ptrtype>>>> M_Fqm_block, M_Lqm_block;

}; // class CRBModelBlock

template <typename SpaceType>
struct InitEnergtMatrixByBlock
{
    typedef std::shared_ptr<SpaceType> space_ptrtype;
    typedef typename Backend<double>::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    explicit InitEnergtMatrixByBlock( space_ptrtype Xh ) :
        m_Xh( Xh )
        {
            m_mats.resize(SpaceType::nSpaces);
        }

    template <typename T>
    void operator()( T const& t ) const
        {
            auto Xh = m_Xh->template functionSpace<T::value>();
            auto u = Xh->element();
            m_mats[T::value] = backend()->newMatrix( _test=Xh, _trial= Xh );
            form2( _test=Xh, _trial=Xh, _matrix=m_mats[T::value] )
                = integrate( _range = elements(Xh->mesh()), _expr = inner(id(u),idt(u)) );
            m_mats[T::value]->close();
        }

    std::vector<sparse_matrix_ptrtype> matrices()
        { return m_mats; }

private:
    space_ptrtype m_Xh;
    mutable std::vector<sparse_matrix_ptrtype> m_mats;
}; // struct InitEnergtMatrixByBlock


template <typename ModelType>
void
CRBModelBlock<ModelType>::initBlocks()
{
    M_backend_l2_vec.resize(n_block);
    for ( int i=0; i<n_block; i++ )
        M_backend_l2_vec[i] = backend( _name="backend-Xh"+std::to_string(i) );

    if ( M_inner_product_matrix_vec.size()==0 )
        M_inner_product_matrix_vec.resize(n_block);

    auto M = this->energyMatrix();
    if ( M )
    {
        Feel::cout << "Using energy matrix as inner product matrix for block CRB !\n";
        for ( int i=0; i<n_block; i++ )
            M_inner_product_matrix_vec[i] = M->createSubMatrix( M->mapRow().dofIdToContainerId(i),
                                                                M->mapCol().dofIdToContainerId(i) );
    }
    else
    {
        InitEnergtMatrixByBlock<space_type> b( this->functionSpace() );
        rangespace_type range;
        boost::fusion::for_each( range, b );
        M_inner_product_matrix_vec = b.matrices();
    }
}


template <typename ModelType>
void
CRBModelBlock<ModelType>::initBlockMatrix()
{
    if ( !M_block_initialized )
    {
        M_block_initialized=true;
        M_Aqm_block.resize( n_block );

        for ( size_type r=0; r<n_block; r++ )
        {
            M_Aqm_block[r].resize( n_block );
            for ( size_type c=0; c<n_block; c++ )
            {
                M_Aqm_block[r][c].resize( this->Qa() );
                for ( size_type q=0; q<this->Qa(); q++ )
                {
                    int mMax = this->mMaxA(q);
                    M_Aqm_block[r][c][q].resize( mMax );
                    for ( size_type m=0; m<mMax; m++ )
                    {
                        auto A = this->M_Aqm[q][m];
                        auto const& i_row = A->mapRow().dofIdToContainerId( r );
                        auto const& i_col = A->mapCol().dofIdToContainerId( c );
                        M_Aqm_block[r][c][q][m] = A->createSubMatrix( i_row, i_col );
                    } // m loop
                } // q loop
            } // c loop
        } // r loop

        int n_output = this->M_Fqm.size();
        M_Fqm_block.resize( n_output );
        for ( size_type l=0; l<n_output; l++ )
        {
            M_Fqm_block[l].resize( n_block );
            for ( size_type r=0; r<n_block; r++ )
            {
                int qMax = this->Ql(l);
                M_Fqm_block[l][r].resize( qMax );
                for ( size_type q=0; q<qMax; q++ )
                {
                    int mMax = this->mMaxF( l, q );
                    M_Fqm_block[l][r][q].resize( mMax );
                    for ( size_type m=0; m<mMax; m++ )
                    {
                        auto V = this->M_Fqm[l][q][m];
                        auto const& i_row = V->map().dofIdToContainerId(r);
                        M_Fqm_block[l][r][q][m] = V->createSubVector( i_row );
                    } // m loop
                } // q loop
            } // block loop
        } // output loop
    } // if !M_block_initialized
}

template <typename ModelType>
void
CRBModelBlock<ModelType>::clearBlockMatrix()
{
    M_Aqm_block.clear();
    M_Fqm_block.clear();
    M_block_initialized = false;
}


} // namespace Feel

#endif
