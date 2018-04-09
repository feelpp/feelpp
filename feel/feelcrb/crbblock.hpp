#ifndef CRBBLOCK_H
#define CRBBLOCK_H

#include <feel/feelcrb/crb.hpp>

namespace Feel
{

po::options_description crbBlockOptions( int const& n_block=3 );

template <typename TruthModelType>
class CRBBlock :
        public CRB<TruthModelType>
{
    typedef CRB<TruthModelType> super;
    typedef CRBBlock<TruthModelType> self_type;


public :
    //@{ // Truth Model
    typedef TruthModelType model_type;
    typedef boost::shared_ptr<model_type> truth_model_ptrtype;
    typedef typename model_type::mesh_type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //@}

    //@{ /// Parameter Space
    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;
    typedef typename parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef typename parameterspace_type::sampling_type sampling_type;
    typedef typename parameterspace_type::sampling_ptrtype sampling_ptrtype;
    //@}

    //@{ /// Function Space and Elements
    typedef typename model_type::space_type space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;
    //@}

    //@{ Backend and Matrix
    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::beta_vector_type beta_vector_type;
    //@}

    //@{ /// Eigen Objects
    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;
    typedef std::vector< std::vector< std::vector< std::vector< matrixN_type >>>> blockmatrixN_type;
    typedef std::vector< std::vector< std::vector< vectorN_type >>> blockvectorN_type;
    //@}

    //@{ /// Exporter
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    //@}

    //@{ /// Database
    typedef CRBElementsDB<model_type> crb_elements_db_type;
    typedef boost::shared_ptr<crb_elements_db_type> crb_elements_db_ptrtype;
    //@}

    typedef typename super::max_error_type max_error_type;
    typedef typename super::error_estimation_type error_estimation_type;
    typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant

    static const int n_block = space_type::nSpaces;
    typedef typename mpl::range_c< int, 0, n_block > rangespace_type;


protected:
    //! Default Constructor
    CRBBlock( std::string const& name = "defaultname_crb",
              crb::stage stage = crb::stage::online,
              WorldComm const& worldComm = Environment::worldComm() ) :
        CRBBlock( name, boost::make_shared<model_type>(stage), stage )
        {}

    //! constructor from command line options
    CRBBlock( std::string const& name, truth_model_ptrtype const & model,
              crb::stage stage = crb::stage::online, std::string const& prefixExt = "" ) :
        super( name, model, stage, prefixExt ),
        M_subN( n_block, 0 )
        {
            for ( int i=0; i<n_block; i++ )
                M_orthonormalize.push_back( boption("crb.block.orthonormalize"+std::to_string(i)) );
        }


public:
    //! \return the expansion of the RB solution \p u, using the block reduced basis
    element_type expansion( vectorN_type const& u, int N = -1,  bool dual=false ) const override;

    void orthonormalizeBasis( int number_of_added_elements ) override;
    bool orthonormalizeSpace( int const& n_space ) const
        {
            CHECK( n_space<n_block ) <<"Invalid space number\n";
            return M_orthonormalize[n_space];
        }

    //! \return the number of basis vectors in space \p n_space
    int& subN( int const& n_space ) { return M_subN[n_space]; }

    std::vector< std::vector< matrixN_type >>& blockAqm( int const& row, int const& col, bool dual=false )
        { return dual ? M_blockAqm_du[row][col]:M_blockAqm_pr[row][col]; }
    std::vector< std::vector< vectorN_type >>& blockFqm( int const& row, bool dual= false )
        { return dual ? M_blockFqm_du[row]:M_blockFqm_pr[row]; }
    std::vector< std::vector< vectorN_type >>& blockLqm( int const& row, bool dual= false )
        { return dual ? M_blockLqm_du[row]:M_blockLqm_pr[row]; }
    std::vector< std::vector< matrixN_type >>& blockAqmPrDu( int const& row, int const& col )
        { return M_blockAqm_pr_du[row][col]; }
    std::vector< std::vector< matrixN_type >>& blockTriqm( int const& row, int const& col )
        {return M_blockTriqm_pr[row][col]; }

    std::vector< std::vector< matrixN_type >> blockAqm( int const& row, int const& col, bool dual=false ) const
        { return dual ? M_blockAqm_du[row][col]:M_blockAqm_pr[row][col]; }
    std::vector< std::vector< vectorN_type >> blockFqm( int const& row, bool dual= false ) const
        { return dual ? M_blockFqm_du[row]:M_blockFqm_pr[row]; }
    std::vector< std::vector< vectorN_type >> blockLqm( int const& row, bool dual= false ) const
        { return dual ? M_blockLqm_du[row]:M_blockLqm_pr[row]; }
    std::vector< std::vector< matrixN_type >> blockAqmPrDu( int const& row, int const& col ) const
        { return M_blockAqm_pr_du[row][col]; }
    std::vector< std::vector< matrixN_type >> blockTriqm( int const& row, int const& col ) const
        {return M_blockTriqm_pr[row][col]; }


    //! save the CRB SP database
    void saveDB() override;
    //! load the CRB SP Database
    bool loadDB() override;

protected:
    void exportBasisFunctions() override;

    void saveRB() override
        { this->M_elements_database.saveDB(); }

    //! add the primal and dual solution \p U and \p Udu to the reduced spaces
    void addBasis( element_type& U, element_type& Udu, parameter_type& mu ) override;

    void buildRbMatrix( int number_of_added_elements, parameter_type& mu, element_ptrtype dual_initial_field ) override;
    void updateAffineDecompositionSize() override;
    void initBlockMatrix();
    virtual void initRezMatrix() {}

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & __ar, const unsigned int __version );

protected:
    std::vector<int> M_subN;
    std::vector<bool> M_orthonormalize;

    blockmatrixN_type M_blockTriqm_pr;
    blockmatrixN_type M_blockAqm_pr, M_blockAqm_du, M_blockAqm_pr_du;
    blockvectorN_type M_blockFqm_pr, M_blockFqm_du, M_blockLqm_pr, M_blockLqm_du;
}; //class CRBBlock


template <typename CRBType>
struct AddBasisByBlock
{
    typedef typename CRBType::element_type element_type;
    typedef typename CRBType::parameter_type parameter_type;

    AddBasisByBlock( element_type& U, element_type& Udu, parameter_type& mu,
                     CRBType* crb ) :
        m_U( U ),
        m_Udu( Udu ),
        m_mu( mu ),
        m_crb( crb )
        {}

    template <typename T>
    void operator()( T const& t ) const
        {
            auto model = m_crb->model();
            auto u = m_U.template elementPtr<T::value>();
            auto udu = m_Udu.template elementPtr<T::value>();

            auto XN = model->rBFunctionSpace()->template rbFunctionSpace<T::value>();
            tic();
            XN->addPrimalBasisElement( u );
            XN->addDualBasisElement( udu );
            m_crb->subN( T::value )++;
            toc("Add Basis Function in space "+std::to_string(T::value) );

            if ( m_crb->model()->addSupremizerInSpace(T::value) )
            {
                tic();
                auto Us = model->supremizer( m_mu, m_U, T::value );
                auto us = Us.template elementPtr<T::value>();
                XN->addPrimalBasisElement( us );
                Us = model->supremizer( m_mu, m_Udu, T::value );
                us = Us.template elementPtr<T::value>();
                XN->addDualBasisElement( us );
                m_crb->subN( T::value )++;
                toc("Add supremizer in space "+std::to_string(T::value) );
            }
        }

private :
    element_type& m_U, m_Udu;
    parameter_type& m_mu;
    CRBType* m_crb;
}; // AddBasisByBlock

template<typename TruthModelType>
void
CRBBlock<TruthModelType>::addBasis( element_type& U, element_type& Udu, parameter_type& mu )
{
    LOG(INFO) << "CRBBlock addBasis begin\n";
    Feel::AddBasisByBlock<self_type> b( U, Udu, mu, this );
    rangespace_type range;
    boost::fusion::for_each( range, b );
     LOG(INFO) << "CRBBlock addBasis end\n";
} // addBasis


template <typename CRBType>
struct OrthonormalizeBasisByBlock
{
    OrthonormalizeBasisByBlock( int& Nm, CRBType* crb ) :
        m_Nm( Nm ),
        m_crb( crb )
        {}

    template <typename T>
    void operator()( T const& t ) const
        {
            if ( m_crb->orthonormalizeSpace(T::value) )
            {
                tic();
                int n_added = m_Nm;
                auto XN = m_crb->model()->rBFunctionSpace()->template rbFunctionSpace<T::value>();
                auto wn = XN->primalRB();
                int N = wn.size();

                if ( m_crb->model()->addSupremizerInSpace(T::value) )
                    n_added *= 2;
                LOG(INFO) <<"CRBBlock orthonomalization begin for space "<<T::value <<", with "<< N
                          <<" basis vectors and "<< n_added <<" new vectors\n";

                double tol = doption(_name="crb.orthonormality-tol");
                int maxit = ioption(_name="crb.orthonormality-max-iter");
                double norm = tol+1;
                int iter=0;
                double old = 0;

                while( norm >=tol && iter < maxit )
                {
                    Feel::cout << "  -- orthonormalization (Gram-Schmidt)\n";
                    for ( size_type i =N-n_added; i < N; ++i )
                    {
                        auto & wni = unwrap_ptr( wn[i] );
                        for ( size_type j = 0; j < i; ++j )
                        {
                            auto & wnj = unwrap_ptr( wn[j] );
                            double __rij_pr = m_crb->model()->scalarProduct(  wni, wnj, T::value );
                            wni.add( -__rij_pr, wnj );
                        }
                    }

                    // normalize
                    for ( size_type i =N-n_added; i < N; ++i )
                    {
                        auto & wni = unwrap_ptr( wn[i] );
                        double __rii_pr = math::sqrt( m_crb->model()->scalarProduct(  wni, wni, T::value ) );
                        wni.scale( 1./__rii_pr );
                    }

                    norm = this->checkOrthonormality( wn, T::value );
                    //if the norm doesn't change
                    if( math::abs(old-norm) < tol )
                        norm=0;
                    old=norm;

                    iter++;
                }
                XN->updatePrimalBasisForUse();
                LOG(INFO) <<"CRBBlock orthonomalization end\n";
                toc( "RB Space Orthonormalization#" + std::to_string(T::value) );
            }
        }

    template <typename WNType>
    double checkOrthonormality ( const WNType& wn, int const& n_space ) const
        {
            int N = wn.size();
            CHECK( N!=0 )<< "[CRB::checkOrthonormality] ERROR : size of wn is zero\n";

            typename CRBType::matrixN_type A, I;
            A.setZero( N, N );
            I.setIdentity( N, N );

            for ( int i = 0; i < N; ++i )
                for ( int j = 0; j < N; ++j )
                    A( i, j ) = m_crb->model()->scalarProduct(  wn[i], wn[j], n_space );

            A -= I;
            LOG( INFO ) << "    o check : " << A.norm() << " (should be 0)";

            return A.norm();
        } //checkOrthonormality()


public:
    int m_Nm;
    CRBType* m_crb;
}; // OrthonormalizeBasisByBlock

template<typename TruthModelType>
void
CRBBlock<TruthModelType>::orthonormalizeBasis( int number_of_added_elements )
{
    LOG(INFO) << "CRBBlock orthonormalize basis begin\n";

    Feel::OrthonormalizeBasisByBlock<self_type> b( number_of_added_elements, this );
    rangespace_type range;
    boost::fusion::for_each( range, b );

    LOG(INFO) << "CRBBlock orthonormalize basis end\n";
}


template <typename CRBType>
struct ExpansionByBlock
{
    typedef typename CRBType::element_type element_type;

    ExpansionByBlock( const CRBType* crb, typename CRBType::vectorN_type const& coeff, int N, bool dual ) :
        m_crb( crb ),
        m_coeff( coeff ),
        m_N( N ),
        m_start( 0 ),
        m_dual( dual )
        {
            m_U = m_crb->model()->functionSpace()->element();
        }

    template <typename T>
    void operator()( T const& t ) const
        {
            auto XN = m_crb->model()->rBFunctionSpace()->template rbFunctionSpace<T::value>();
            auto WN = m_dual ? XN->dualRB() : XN->primalRB();
            int Nwn = m_N>0 ? m_N:m_crb->dimension();
            if ( m_crb->model()->addSupremizerInSpace(T::value) )
                Nwn *= 2;
            CHECK( Nwn<=WN.size() ) << "invalid expansion size\n";

            CHECK( m_start+Nwn<=m_coeff.size() ) << "invalide expansion size\n";
            auto  coeff = m_coeff.segment( m_start, Nwn );
            auto u = m_U.template element<T::value>();
            u = Feel::expansion( WN, coeff, Nwn ).container();
            m_start += Nwn;
        }

    typename CRBType::element_type U() const { return m_U; }

private:
    const CRBType* m_crb;
    const typename CRBType::vectorN_type& m_coeff;
    int m_N;
    mutable int m_start;
    bool m_dual;
    element_type m_U;
}; // ExpansionByBlock

template<typename TruthModelType>
typename CRBBlock<TruthModelType>::element_type
CRBBlock<TruthModelType>::expansion( vectorN_type const& u, int N, bool dual ) const
{
    LOG(INFO) << "Expansion by block begin\n";

    Feel::ExpansionByBlock<self_type> b( this, u, N, dual );
    rangespace_type range;
    boost::fusion::for_each( range, b );

    LOG(INFO) << "Expansion by block end\n";
    return b.U();
}


template <typename CRBType>
struct BuildRbMatrixByRow
{
    typedef typename CRBType::rangespace_type range_type;

    BuildRbMatrixByRow( CRBType* crb, int const& n_added ) :
        m_crb( crb ),
        m_n_added( n_added )
        {}

    template <typename R>
    void operator()( R const& t ) const
        {
            BuildRbMatrixByCol<R> col( m_crb, m_n_added );
            range_type range;
            boost::fusion::for_each( range, col );
        }

    template <typename R>
    struct BuildRbMatrixByCol
    {
        typedef typename CRBType::element_type element_type;

        BuildRbMatrixByCol( CRBType* crb, int const& n_added ) :
            m_crb( crb ),
            m_n_added( n_added )
            {}

        template <typename C>
        void operator()( C const& t ) const
            {
                LOG(INFO) << "CRBBlock assemble RB matrices block "<<R::value<<C::value <<std::endl;

                auto model = m_crb->model();

                auto XNr = model->rBFunctionSpace()->template rbFunctionSpace<R::value>();
                auto XNc = model->rBFunctionSpace()->template rbFunctionSpace<C::value>();

                element_type Ur = model->functionSpace()->element();
                element_type Uc = model->functionSpace()->element();

                auto ur = Ur.template element<R::value>();
                auto uc = Uc.template element<C::value>();

                int Nr = m_crb->subN(R::value);
                int Nc = m_crb->subN(C::value);

                int n_upr = m_n_added;
                int n_upc = m_n_added;
                if ( model->addSupremizerInSpace(R::value) )
                    n_upr *= 2;
                if ( model->addSupremizerInSpace(C::value) )
                    n_upc *= 2;
                if( ioption(_name="ser.rb-frequency")!=0 && !m_crb->rebuild() )
                {
                    n_upr = Nr;
                    n_upc = Nc;
                }


                for ( size_type q=0; q<model->Qa(); q++ )
                {
                    // resize the block RC
                    for ( size_type m=0; m<model->mMaxA(q); m++ )
                    {
                        m_crb->blockAqm( R::value, C::value, false )[q][m].conservativeResize( Nr, Nc );
                        m_crb->blockAqm( R::value, C::value, true )[q][m].conservativeResize( Nr, Nc );
                        m_crb->blockAqmPrDu( R::value, C::value )[q][m].conservativeResize( Nr, Nc );

                        // update last rows of block RC
                        for ( size_type i= Nr-n_upr; i<Nr; i++ )
                        {
                            for ( size_type j=0; j<Nc; j++ )
                            {
                                Ur.zero();
                                ur = XNr->primalBasisElement(i);
                                Uc.zero();
                                uc = XNc->primalBasisElement(j);
                                m_crb->blockAqm( R::value, C::value, false )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur );
                                ur = XNr->dualBasisElement(i);
                                uc = XNc->dualBasisElement(j);
                                m_crb->blockAqm( R::value, C::value, true )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur, true );
                                ur = XNr->dualBasisElement(i);
                                uc = XNc->primalBasisElement(j);
                                m_crb->blockAqmPrDu( R::value, C::value )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur );
                            } // loop on j
                        } // loop on i

                        // update last col of block RC
                        for ( size_type i=0; i<Nr; i++ )
                        {
                            for ( size_type j=Nc-n_upc; j<Nc; j++ )
                            {
                                ur = XNr->primalBasisElement(i);
                                uc = XNc->primalBasisElement(j);
                                m_crb->blockAqm( R::value, C::value, false )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur );
                                ur = XNr->dualBasisElement(i);
                                uc = XNc->dualBasisElement(j);
                                m_crb->blockAqm( R::value, C::value, true )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur, true );
                                ur = XNr->dualBasisElement(i);
                                uc = XNc->primalBasisElement(j);
                                m_crb->blockAqmPrDu( R::value, C::value )[q][m](i,j)
                                    = model->Aqm( q, m, Uc, Ur );
                            } // loop on j
                        } // loop on i
                    } // loop on m
                } //loop on q

                // update rhs blocks if c==0
                if ( C::value==0 )
                {
                    for ( size_type q=0; q<model->Ql(0); q++ )
                    {
                        for ( size_type m=0; m<model->mMaxF(0, q); m++ )
                        {
                            m_crb->blockFqm( R::value, false )[q][m].conservativeResize(Nr);
                            m_crb->blockFqm( R::value, true )[q][m].conservativeResize(Nr);

                            for ( size_type l=Nr-n_upr; l<Nr; l++ )
                            {
                                Ur.zero();
                                ur=XNr->primalBasisElement(l);
                                m_crb->blockFqm( R::value, false )[q][m](l) = model->Fqm( 0, q, m, Ur );
                                ur=XNr->dualBasisElement(l);
                                m_crb->blockFqm( R::value, true )[q][m](l) = model->Fqm( 0, q, m, Ur );
                            } // loop on l
                        } // loop on m
                    } //loop on q

                    int oi = m_crb->outputIndex();
                    for ( size_type q=0; q<model->Ql(oi); q++ )
                    {
                        for ( size_type m=0; m<model->mMaxF(oi, q); m++ )
                        {
                            m_crb->blockLqm( R::value, false )[q][m].conservativeResize(Nr);
                            m_crb->blockLqm( R::value, true )[q][m].conservativeResize(Nr);

                            for ( size_type l=Nr-n_upr; l<Nr; l++ )
                            {
                                Ur.zero();
                                ur=XNr->primalBasisElement(l);
                                m_crb->blockLqm( R::value, false )[q][m](l) = model->Fqm( oi, q, m, Ur );
                                ur=XNr->dualBasisElement(l);
                                m_crb->blockLqm( R::value, true )[q][m](l) = model->Fqm( oi, q, m, Ur );
                            } // loop on l
                        } // looop on m
                    } //loop on q
                } // if c==0
            }

    private:
        CRBType* m_crb;
        int m_n_added;
    };

private:
    CRBType* m_crb;
    int m_n_added;
}; //struct BuildRbMatrixByRow


template<typename TruthModelType>
void
CRBBlock<TruthModelType>::buildRbMatrix( int number_of_added_elements, parameter_type& mu,
                                         element_ptrtype dual_initial_field )
{
    tic();
    LOG(INFO) <<"CRBBlock Build rb matrices begin\n";

    Feel::BuildRbMatrixByRow<self_type> row( this, number_of_added_elements );
    rangespace_type range;
    boost::fusion::for_each( range, row );

    LOG(INFO) <<"CRBBlock Build rb matrices end\n";
    toc("Reduced Matrices Built");
}// buildRbMatrix


template<typename TruthModelType>
void
CRBBlock<TruthModelType>::initBlockMatrix()
{
    M_blockAqm_pr.resize(n_block);
    M_blockAqm_du.resize(n_block);
    M_blockAqm_pr_du.resize(n_block);
    for( int i = 0; i < n_block; i++ )
    {
        M_blockAqm_pr[i].resize(n_block);
        M_blockAqm_du[i].resize(n_block);
        M_blockAqm_pr_du[i].resize(n_block);
    }

    M_blockFqm_pr.resize(n_block);
    M_blockLqm_pr.resize(n_block);
    M_blockFqm_du.resize(n_block);
    M_blockLqm_du.resize(n_block);

    int q_tri = this->M_model->QTri();
    if ( q_tri )
    {
        M_blockTriqm_pr.resize( n_block );
        for( int i = 0; i < n_block; i++ )
            M_blockTriqm_pr[i].resize(n_block);

    }
}

template<typename TruthModelType>
void
CRBBlock<TruthModelType>::updateAffineDecompositionSize()
{
    if( this->M_rebuild || this->M_N == 0 )
        initBlockMatrix();
    int output_index = this->M_output_index;

    for ( int n=0; n<n_block; n++ )
        subN(n) = this->M_model->addSupremizerInSpace(n) ? 2*this->M_N:this->M_N;

    for (int r=0; r<n_block; r++ )
    {
        for ( int c=0; c<n_block; c++ )
        {
            M_blockAqm_pr[r][c].resize( this->M_model->Qa() );
            M_blockAqm_du[r][c].resize( this->M_model->Qa() );
            M_blockAqm_pr_du[r][c].resize( this->M_model->Qa() );
            for ( int q=0; q<M_blockAqm_pr[r][c].size(); q++ )
            {
                M_blockAqm_pr[r][c][q].resize( this->M_model->mMaxA(q) );
                M_blockAqm_du[r][c][q].resize( this->M_model->mMaxA(q) );
                M_blockAqm_pr_du[r][c][q].resize( this->M_model->mMaxA(q) );
            }
        }

        M_blockFqm_pr[r].resize( this->M_model->Ql(0) );
        M_blockFqm_du[r].resize( this->M_model->Ql(0) );
        for ( int q=0; q<M_blockFqm_pr[r].size(); q++ )
        {
            M_blockFqm_pr[r][q].resize( this->M_model->mMaxF( 0, q) );
            M_blockFqm_du[r][q].resize( this->M_model->mMaxF( 0, q) );
        }

        M_blockLqm_pr[r].resize( this->M_model->Ql( output_index ) );
        M_blockLqm_du[r].resize( this->M_model->Ql( output_index ) );
        for ( int q=0; q<M_blockLqm_pr[r].size(); q++ )
        {
            M_blockLqm_pr[r][q].resize( this->M_model->mMaxF( output_index, q ) );
            M_blockLqm_du[r][q].resize( this->M_model->mMaxF( output_index, q ) );
        }
    }

    int q_tri = this->M_model->QTri();
    if ( q_tri )
    {
        for (int r=0; r<n_block; r++ )
            for ( int c=0; c<n_block; c++ )
                M_blockTriqm_pr[r][c].resize(q_tri);
    }

    if ( this->M_error_type == CRB_RESIDUAL || this->M_error_type == CRB_RESIDUAL_SCM )
        initRezMatrix();

} // updateAffinedecompositionsize()



template <typename CRBType>
struct ExportBasisFunctionsByBlock
{
    typedef typename CRBType::mesh_type mesh_type;

    ExportBasisFunctionsByBlock( CRBType* crb ):
        m_crb( crb )
        {
            m_e = exporter( _mesh=m_crb->model()->functionSpace()->mesh(), _name="basis-functions" );
        }

    template <typename T>
    void operator()( T const& t ) const
        {
            auto model = m_crb->model();
            auto XN = model->rBFunctionSpace()->template rbFunctionSpace<T::value>();
            auto u_vec = XN->primalRB();

            CHECK( u_vec.size()>0 )<< "CRBBlock error : no element to export in space "<<T::value <<std::endl;

            for ( int index=0; index<u_vec.size(); index++ )
            {
                int n = m_crb->model()->addSupremizerInSpace(T::value) ? index/2:index;
                auto mu = m_crb->wnmu()->at( n );

                std::string basis_name = ( boost::format("u%1%_pr_%2%_param")%T::value %index ).str();
                std::string name = basis_name + mu.toString();
                m_e->add( name, unwrap_ptr( u_vec[index] ) );
            }
        }

    void save() { m_e->save(); }

private:
    CRBType* m_crb;
    boost::shared_ptr<Exporter<mesh_type>> m_e;

}; // struct ExportBasisFunctionsByBlock

template<typename TruthModelType>
void
CRBBlock<TruthModelType>::exportBasisFunctions()
{
    tic();
    LOG(INFO) << "CRBBLock exportbasis functions begin\n";
    ExportBasisFunctionsByBlock<self_type> b( this );
    rangespace_type range;
    boost::fusion::for_each( range, b );
    b.save();
    LOG(INFO) << "CRBBLock exportbasis functions end\n";
    toc("CRBBlock Export basis functions");
}




template<typename TruthModelType>
template<class Archive>
void
CRBBlock<TruthModelType>::serialize(Archive & ar, const unsigned int __version )
{
    ar & boost::serialization::base_object<super>( *this );

    ar & BOOST_SERIALIZATION_NVP( M_subN );

    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_pr );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockFqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockLqm_du );
    ar & BOOST_SERIALIZATION_NVP( M_blockAqm_pr_du );
}

template<typename TruthModelType>
void
CRBBlock<TruthModelType>::saveDB()
{
    fs::ofstream ofs( this->dbLocalPath() / this->dbFilename() );

    if ( ofs )
    {
        // boost::archive::text_oarchive oa( ofs );
        boost::archive::binary_oarchive oa( ofs );
        // write class instance to archive
        oa << *this;
        // archive and stream closed when destructors are called
    }
    this->saveJson();
} //saveDB()

template<typename TruthModelType>
bool
CRBBlock<TruthModelType>::loadDB()
{
    if( this->isDBLoaded() )
        return true;

    fs::path db = this->lookForDB();

    if ( db.empty() )
        return false;

    if ( !fs::exists( db ) )
        return false;

    fs::ifstream ifs( db );

    if ( ifs )
    {
        // boost::archive::text_iarchive ia( ifs );
        boost::archive::binary_iarchive ia( ifs );
        // write class instance to archive
        ia >> *this;
        this->setIsLoaded( true );
        // archive and stream closed when destructors are called
        return true;
    }

    return false;
} // loadDB()







} //namespace Feel

#endif
