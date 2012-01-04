/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2008-01-03

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file matrixblock.cpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2011-06-10
 */

//#include <feel/feelalg/matrixblock.hpp>



namespace Feel
{


template <int NR, int NC,typename T>
MatrixBlock<NR,NC,T>::MatrixBlock( vf::Blocks<NR,NC,matrix_ptrtype> const & blockSet,
                                   backend_type &backend,
                                   bool _doAssemble)
    :
    super(),
    //M_v( v ),
    M_mat()
{

    auto v = blockSet.getSetOfBlocks();

    //std::cout << "[MatrixBlock::MatrixBlock] compute size1 size2" << std::endl;
    size_type _size2 = 0;
    for (uint i=0;i<NC;++i)
        _size2 += v[i]->size2();

    size_type _size1 =0;
    for (uint i=0;i<NR;++i)
        _size1 += v[i*NC]->size1();

    //std::cout << "[MatrixBlock::MatrixBlock] build graph" << std::endl;
    graph_ptrtype graph(new graph_type(0,0,_size1-1,0,_size2-1));//( new graph_type( ) );
    size_type start_i=0;
    size_type start_j=0;
    for (uint i=0;i<NR;++i)
        {
            start_j=0;
            for (uint j=0;j<NC;++j)
                {
                    v[i*NC+j]->close();
                    this->mergeBlockGraph(graph,v[i*NC+j],start_i,start_j);

                    start_j += v[i*NC+j]->size2();
                }
            start_i += v[i*NC]->size1();
        }

    graph->close();

    // std::cout << "[MatrixBlock::MatrixBlock] build Matrix" << std::endl;
    M_mat = backend.newMatrix(_size1,_size2,_size1,_size2,graph);
    M_mat->zero();
    //M_mat->graph()->showMe();

    if (_doAssemble)
        {
            start_i=0;
            start_j=0;
            for (uint i=0;i<NR;++i)
                {
                    start_j=0;
                    for (uint j=0;j<NC;++j)
                        {
                            this->updateBlockMat(v[i*NC+j],start_i,start_j);

                            start_j += v[i*NC+j]->size2();
                        }
                    start_i += v[i*NC]->size1();
                }
        }


    //std::cout << "[MatrixBlock::MatrixBlock] build FieldSplit index" << std::endl;
    // index container for field split preconditioner
    std::vector < std::vector<int> > indexSplit(NR);
    uint16_type startIS = 0;
    for (uint i=0;i<NR;++i)
        {
            auto Loc_nDof = v[i*NC]->size1();
            indexSplit[i].resize(Loc_nDof);
            for (uint l = 0; l< Loc_nDof ; ++l )
                {
                    indexSplit[i][l] = startIS + l;
                }
            startIS += Loc_nDof;
        }
    // update
    M_mat->setIndexSplit(indexSplit);

}


template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::mergeBlockGraph(graph_ptrtype & globGraph,
                                      matrix_ptrtype m,
                                      size_type start_i, size_type start_j)
{

    //FEEL_ASSERT(m->hasGraph()).error("sub matrix doesn t have a graph");

    auto g = m->graph();

    auto it = g->begin();
    auto en = g->end();
    for( ; it != en; ++it )
        {
            int theglobalrow = start_i + it->first;
            typename graph_type::row_type & row = globGraph->row(theglobalrow);

            //globGraph->row(theglobalrow).template get<0>() = it->second.template get<0>();//rank
            row.template get<0>() = it->second.template get<0>();//rank

            int thelocalrow = start_i + it->second.template get<1>();
            //globGraph->row(theglobalrow).template get<1>() = thelocalrow;
            row.template get<1>() = thelocalrow;

            auto nbDof = it->second.template get<2>().size();
            if (nbDof>0)
                {
                    // Get the row of the sparsity pattern
#if 0
                    std::vector<size_type> ivec(  it->second.template get<2>().begin(),  it->second.template get<2>().end() );
                    std::for_each( ivec.begin(), ivec.end(), boost::phoenix::arg_names::arg1 += start_j );
#else

                    std::vector<size_type> ivec(nbDof);
                    auto itDof=it->second.template get<2>().begin();
                    for ( int i=0;
                          i<nbDof;
                          ++i,++itDof) ivec[i]=*itDof+start_j;
#endif
                    //std::set<size_type> iout( ivec.size()+ M_graph->row(theglobalrow).template get<2>().size() );
                    //std::set<size_type> iout( ivec.begin(), ivec.end() );
                    //iout.insert( globGraph->row(theglobalrow).template get<2>().begin(),
                    //             globGraph->row(theglobalrow).template get<2>().end() );
                    //globGraph->row(theglobalrow).template get<2>() = iout;
                    //globGraph->row(theglobalrow).template get<2>().insert(ivec.begin(), ivec.end());
                    row.template get<2>().insert(ivec.begin(), ivec.end());
                }

        }


}



template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::init (const size_type m,
                            const size_type n,
                            const size_type m_l,
                            const size_type n_l,
                            const size_type nnz,
                            const size_type /*noz*/)
{
    if ((m==0) || (n==0))
        return;

    {
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        this->setInitialized( true );
    }

    M_mat->init(m,n,m_l,n_l,nnz);

}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::init (const size_type m,
                           const size_type n,
                           const size_type m_l,
                           const size_type n_l,
                           graph_ptrtype const& graph )
{
    this->setGraph( graph );

    {
        // Clear initialized matrices
        if (this->isInitialized())
            this->clear();

        this->setInitialized(  true );
    }

    M_mat->init(m,n,m_l,n_l,graph);
}


template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::zero ()
{
    M_mat->zero();
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::zero ( size_type /*start1*/, size_type /*stop1*/, size_type /*start2*/, size_type /*stop2*/ )
{
    M_mat->zero();
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::clear ()
{
    M_mat->clear();
}

template <int NR, int NC,typename T>
inline
void
MatrixBlock<NR,NC,T>::close () const
{
    M_mat->close();
}

template <int NR, int NC,typename T>
inline
size_type
MatrixBlock<NR,NC,T>::size1 () const
{
    return M_mat->size1();
}

template <int NR, int NC,typename T>
inline
size_type
MatrixBlock<NR,NC,T>::size2 () const
{
    return M_mat->size2();
}

template <int NR, int NC,typename T>
inline
size_type
MatrixBlock<NR,NC,T>::rowStart () const
{
    return M_mat->rowStart();
}

template <int NR, int NC,typename T>
inline
size_type
MatrixBlock<NR,NC,T>::rowStop () const
{
    return M_mat->rowStop();
}

template <int NR, int NC,typename T>
inline
void
MatrixBlock<NR,NC,T>::set (const size_type i,
                           const size_type j,
                           const value_type& value)
{
    M_mat->set(i,j,value);
}

template <int NR, int NC,typename T>
inline
void
MatrixBlock<NR,NC,T>::add (const size_type i,
                           const size_type j,
                           const value_type& value)
{
    M_mat->add(i,j,value);
}

template <int NR, int NC,typename T>
inline
bool
MatrixBlock<NR,NC,T>::closed() const
{
    return M_mat->closed();
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T> ::addMatrix(const ublas::matrix<value_type>& dm,
                                 const std::vector<size_type>& rows,
                                 const std::vector<size_type>& cols)
{
    M_mat->addMatrix(dm,rows,cols);
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::addMatrix ( int* rows, int nrows,
                                  int* cols, int ncols,
                                  value_type* data )
{
    M_mat->addMatrix(rows,nrows,cols,ncols,data);
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::printMatlab (const std::string name) const
{
    M_mat->printMatlab(name);
}

template <int NR, int NC,typename T>
inline
void
MatrixBlock<NR,NC,T>::addMatrix (const value_type a_in, MatrixSparse<value_type> &X_in)
{
    M_mat->addMatrix(a_in,X_in);
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::scale( value_type const a )
{
    M_mat->scale(a);
}

template <int NR, int NC,typename T>
typename MatrixBlock<NR,NC,T>::value_type
MatrixBlock<NR,NC,T>::energy( Vector<value_type> const& __v,
                              Vector<value_type> const& __u,
                              bool _transpose ) const
{
    return M_mat->energy(__v,__u,_transpose);
}

template <int NR, int NC,typename T>
inline
typename MatrixBlock<NR,NC,T>::real_type
MatrixBlock<NR,NC,T>::l1Norm() const
{
    return M_mat->l1Norm();
}

template <int NR, int NC,typename T>
inline
typename MatrixBlock<NR,NC,T>::real_type
MatrixBlock<NR,NC,T>::linftyNorm() const
{
    return M_mat->linftyNorm();
}

template <int NR, int NC,typename T>
inline
typename MatrixBlock<NR,NC,T>::value_type
MatrixBlock<NR,NC,T>::operator () (const size_type i,
                                   const size_type j) const
{
    return (*M_mat)(i,j);
}


template <int NR, int NC,typename T>
MatrixBlock<NR,NC,T> &
MatrixBlock<NR,NC,T>::operator = ( MatrixSparse<value_type> const& M )
{
    *M_mat = M;
    return *this;
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::zeroRows( std::vector<int> const& rows, std::vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context )
{
    M_mat->zeroRows(rows,values,rhs,on_context);
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::transpose( MatrixSparse<value_type>& Mt ) const
{
    M_mat->transpose(Mt);
}

template <int NR, int NC,typename T>
void
MatrixBlock<NR,NC,T>::updateBlockMat(boost::shared_ptr<MatrixSparse<value_type> > m, size_type start_i, size_type start_j)
{
        M_mat->updateBlockMat(m,start_i,start_j);
}


} // Feel


