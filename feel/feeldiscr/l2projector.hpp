#ifndef _L2PROJECTOR_HPP_
#define _L2PROJECTOR_HPP_

#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

template<class DomainSpace, class DualImageSpace>
class L2Projector : public OperatorLinear<DomainSpace, DualImageSpace>
{
    typedef L2Projector<DomainSpace,DualImageSpace> super;

public :

    // -- TYPEDEFS --
  // typedef Operator<DomainSpace, DualImageSpace> super_type;
    typedef OperatorLinear<DomainSpace, DualImageSpace> ol_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::dual_image_space_type  dual_image_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename super::dual_image_space_ptrtype  dual_image_space_ptrtype;
    typedef typename domain_space_type::element_type domain_element_type;
    typedef typename dual_image_space_type::element_type dual_image_element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;


  L2Projector(domain_space_ptrtype     domainSpace,
	       dual_image_space_ptrtype dualImageSpace,
	       backend_ptrtype          backend)
    :
    ol_type(domainSpace, dualImageSpace, backend),
    M_backend(backend),
    M_matrix(M_backend->newMatrix( domainSpace, dualImageSpace ))
  {  initMatrix();  }

  // //used in the case where domainSpace and dualImageSpace are the same
  L2Projector(domain_space_ptrtype     domainSpace,
  	       backend_ptrtype          backend)
    :
    ol_type(domainSpace, domainSpace, backend),
    M_backend(backend),
    M_matrix(M_backend->newMatrix( domainSpace, domainSpace ))
  {  initMatrix();  }

  ~L2Projector() {}

  template<typename RhsExpr>
  domain_element_type
  project( RhsExpr const& rhs_expr )
  {
    domain_element_type de = this->domainSpace()->element();
    
    auto ie = M_backend->newVector(this->dualImageSpace());
    form1(_test=this->dualImageSpace(), _vector=ie, _init=true) =
      integrate(elements(this->domainSpace()->mesh()),
		rhs_expr * id( this->dualImageSpace()->element() ) );

    M_backend->solve(M_matrix, de, ie);

    return de;
  }


  template<typename RhsExpr>
  domain_element_type
  operator()(RhsExpr const& rhs_expr)
  {   return this->project(rhs_expr);   }



  template<typename RhsExpr>
  void
  operator()(domain_element_type& de, RhsExpr const& rhs_expr)
  {   de = this->project(rhs_expr);   }



  domain_element_type 
  operator()(image_element_type const& ie)
  {   
    domain_element_type de = this->domainSpace()->element();
    M_backend->solve(M_matrix, de, ie);
    return de ;
  }



  void 
  operator()(domain_element_type &de, image_element_type const& ie)
  {   M_backend->solve(M_matrix, de, ie);  }



  void 
  apply(domain_element_type& de,
	     image_element_type const& ie)
  {     M_backend->solve(M_matrix, de, ie);    }



  template<typename RhsExpr>
  void
  apply(domain_element_type& de,
	RhsExpr const& rhs_expr)
  {    de=this->project(rhs_expr);  }


#if 0
  //a finir
template<typename TDomainSpace, typename TDualImageSpace>
boost::shared_ptr< L2Projector<TDomainSpace, TDualImageSpace> >
l2Projector( boost::shared_ptr<TDomainSpace> const& domainspace,
	     boost::shared_ptr<TDualImageSpace> const& imagespace,
	     typename L2Projector<TDomainSpace, TDualImageSpace>::backend_ptrtype const& backend )
  {
    typedef L2Projector<TDomainSpace, TDualImageSpace> L2Proj_type;
    boost::shared_ptr<L2Proj_type> l2p( new L2Proj_type(domainspace, imagespace, backend) );
    return l2p;
  }
#endif



private :

  void initMatrix()
  {
    form2 (_trial=this->domainSpace(),
	   _test=this->dualImageSpace(),
	   _matrix=M_matrix) =
      integrate(elements(this->domainSpace()->mesh()),
		trans(idt( this->domainSpace()->element() )) /*trial*/
		*id( this->domainSpace()->element() ) /*test*/
		);
    M_matrix->close();
  }

  backend_ptrtype M_backend;
  matrix_ptrtype M_matrix;

};//L2Projector

} //namespace Feel


#endif
