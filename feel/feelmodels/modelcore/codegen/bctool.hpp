#ifndef __CLDATA_H
#define __CLDATA_H 1

#include <boost/mpl/list_c.hpp>
#include <boost/mpl/push_front.hpp>
#include <boost/mpl/count_if.hpp>
//#include <feel/feelcore/feel.hpp>
#include <feel/feelvf/vf.hpp>

//#include <boost/phoenix.hpp>
//#include <research/life-workspace/Fluid-Structure/cltype.hpp>
//#include <fsi/fsicore/bcmarker.cpp>

#ifndef ADD_MARKER
#define ADD_MARKER(s) (s,#s)
#endif

//#if defined( FSI_USE_DEFAULT_MARKER )
//#include <fsi/fsicore/default.bcmarkers>
//#else
//#include BC_MARKER_NAME
#include "bcmarker.cpp"
//#endif
namespace Feel
{


#define DIRICHLET_TYPE_NAME2(i)                                     \
    BOOST_PP_TUPLE_ELEM(2, 0, BOOST_PP_ARRAY_ELEM(i,MARKER_LIST)) \
    /**/
#define DIRICHLET_TYPE_STR(i)                                       \
    BOOST_PP_TUPLE_ELEM(2, 1, BOOST_PP_ARRAY_ELEM(i,MARKER_LIST)) \


    /////////////////////////////////////////////////////////////////////////

#define CL_INSTANTIATES_FOR_COMP(r, state)                              \
    BOOST_PP_NOT_EQUAL( BOOST_PP_TUPLE_ELEM(2, 0, state),               \
                        BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 1, state))  \
                        )                                               \
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define CL_INSTANTIATES_FOR_INCR(r, state)              \
    ( BOOST_PP_INC(BOOST_PP_TUPLE_ELEM(2, 0, state)),   \
      BOOST_PP_TUPLE_ELEM(2, 1, state) )                \
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define CL_INSTANTIATES_FOR(r,state)                                    \
    DIRICHLET_TYPE_NAME2(BOOST_PP_TUPLE_ELEM(2,0,state)) BOOST_PP_COMMA() \
    /**/
    /*_________________________________________________*/
    /*                                                 */
    /**/
#define CL_MARKER_FOR(r,state)                                          \
    DIRICHLET_TYPE_STR(BOOST_PP_TUPLE_ELEM(2,0,state))  BOOST_PP_COMMA() \
    /*_________________________________________________*/
    /*                                                 */
    /**/



    namespace cl {

        enum cl_markerName_enum
            {
                BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(MARKER_LIST), 2) ),
                              CL_INSTANTIATES_FOR_COMP,
                              CL_INSTANTIATES_FOR_INCR,
                              CL_INSTANTIATES_FOR )
                /**/
                DIRICHLET_TYPE_NAME2(BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(MARKER_LIST), 1))
                /**/
            };

        static const char * cl_marker_string[] =
            {
                /**/
                BOOST_PP_FOR( (0, BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(MARKER_LIST), 2) ),
                              CL_INSTANTIATES_FOR_COMP,
                              CL_INSTANTIATES_FOR_INCR,
                              CL_MARKER_FOR )
                /**/
                DIRICHLET_TYPE_STR(BOOST_PP_SUB(BOOST_PP_ARRAY_SIZE(MARKER_LIST), 1))
                /**/
            };


        enum cl_type_enum
            {
                dirichlet=0,
                dirichlet_vec,
                dirichlet_x,
                dirichlet_y,
                dirichlet_z,
                paroi_mobile,
                neumann_scal,
                neumann_vec,
                slip,
                pressure,
                fbm,
                fbm_dirichlet,
                fbm_dirichlet_vec,
                fbm_moving_wall,
                fluid_outlet,
                robin_vec,
                follower_pressure
            };
#if 0
        static const char * cl_type_string[] =
            {
                "dirichlet",
                "dirichlet_vect",
                /*"dirichlet_x",
                "dirichlet_y",
                 "dirichlet_z",*/
                "paroi_mobile",
                "neumann_scal",
                "neumann_vec",
                "slip",
                "pressure",
                "fbm",
                "fbm_dirichlet",
                "fbm_dirichlet_vec",
                "fbm_moving_wall",
                "fluid_outlet",
                "robin_vec",
                "follower_pressure"
            };
#endif



    } // end namespace cl

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    //Macros user
#if 1
#define ForEachBC(CLBASE,CLT,EXPR)                 \
    for (int i=0;i<CLBASE.nbCL<CLT>();++i)                              \
        {                                                               \
            switch (i)                                                  \
                {                                                       \
                case 0 :  {                                             \
                    auto PhysicalName = CLBASE.getMarkerNameByLocalBcId< CLT,0 >(); \
                    /*auto numExpr = CLBASE.localToGlobal<CLT,0>();*/   \
                    auto Expression = CLBASE.getExprByLocalBcId<CLT,0>(); \
                    EXPR;                                               \
                    break; }                                            \
                case 1 :  {                                             \
                    auto PhysicalName = CLBASE.getMarkerNameByLocalBcId< CLT,1 >(); \
                    /*auto numExpr = CLBASE.localToGlobal<CLT,1>();*/   \
                    auto Expression = CLBASE.getExprByLocalBcId<CLT,1>();\
                    EXPR;                                               \
                    break; }                                            \
                case 2 :  {                                             \
                    auto PhysicalName = CLBASE.getMarkerNameByLocalBcId< CLT,2 >(); \
                    /*auto numExpr = CLBASE.localToGlobal<CLT,2>();*/   \
                    auto Expression = CLBASE.getExprByLocalBcId<CLT,2>();\
                    EXPR;                                               \
                    break;}                                             \
                case 3 :  {                                             \
                    auto PhysicalName = CLBASE.getMarkerNameByLocalBcId< CLT,3 >(); \
                    /*auto numExpr = CLBASE.localToGlobal<CLT,3>();*/   \
                    auto Expression = CLBASE.getExprByLocalBcId<CLT,3>(); \
                    EXPR;                                               \
                    break;}                                             \
                case 4 :  {                                             \
                    auto PhysicalName = CLBASE.getMarkerNameByLocalBcId< CLT,4 >(); \
                    /*auto numExpr = CLBASE.localToGlobal<CLT,4>();*/   \
                    auto Expression = CLBASE.getExprByLocalBcId<CLT,4>();\
                    EXPR;                                               \
                    break;}                                             \
                }                                                       \
        }                                                               \
    /**/
#else

        template< typename TheBcDataBaseType>
        struct ForEachBCImpl
        {
            ForEachBCImpl( TheBcDataBaseType const& bcData, std::function<void(std::string,LambdaExpr1)> func )
                :
                M_bcData( bcData ),
                M_function( func )
                {}

            template<typename TheBcDataType>
            void operator()( TheBcDataType const& _bcDataType ) const
            {
                std::string PhysicalName = M_bcData.template getMarkerName< TheBcDataType::value >();
                LambdaExpr1 hola;
                M_function( PhysicalName,hola/*LambdaExpr1*//*Feel::vf::_e1*/ );
            }
        private :
            TheBcDataBaseType const& M_bcData;
            std::function<void(std::string,LambdaExpr1)> M_function;
        };

#define ForEachBC(CLBASE,CLT,EXPR)                                      \
    {                                                                   \
        std::function<void (int)> myfunc = [&](int kkk) { int oo = CLBASE.nbCL<CLT>(); }; \
        auto const& TheLambdaExpression = vf::vec(Feel::vf::_e1,Feel::vf::_e2,Feel::vf::_e3); \
        std::function<void(std::string,LambdaExpr1)> myfunc2 = [&](std::string PhysicalName,LambdaExpr1 const& mylambda1/**/) \
            {                                                           \
                /*auto Expression = vf::vec(Feel::vf::_e1,Feel::vf::_e2,Feel::vf::_e3);*/ \
                auto Expression = TheLambdaExpression;                  \
                /*auto Expression = Feel::vf::_e1;*/                    \
                EXPR;                                                   \
            };                                                          \
                                                                        \
        auto myconnntaine = CLBASE.globalIndicesByBc<CLT>();            \
        Feel::ForEachBCImpl<decltype(CLBASE)/*,decltype(TheLambdaExpression)*/> myImpl(CLBASE,myfunc2); \
        fusion::for_each( myconnntaine, myImpl );                       \
        /*fusion::for_each( myconnntaine, [&]() { auto PhysicalName = CLBASE.getMarkerName< boost::phoenix::placeholders::arg1 >(); } );*/ \
        double hollla=3.2;                                              \
    }                                                                   \
        /**/
#endif


#define ForEachBC_TEMPLATE(CLBASE,CLT,EXPR)                             \
    for (int i=0;i<CLBASE.template nbCL<CLT>();++i)                     \
        {                                                               \
            switch (i)                                                  \
                {                                                       \
                case 0 :  {                                             \
                    auto PhysicalName = CLBASE.template getMarkerNameByLocalBcId< CLT,0 >(); \
                    auto Expression = CLBASE.template getExprByLocalBcId<CLT,0>(); \
                    EXPR;                                               \
                    break; }                                            \
                case 1 :  {                                             \
                    auto PhysicalName = CLBASE.template getMarkerNameByLocalBcId< CLT,1 >(); \
                    auto Expression = CLBASE.template getExprByLocalBcId<CLT,1>();\
                    EXPR;                                               \
                    break; }                                            \
                case 2 :  {                                             \
                    auto PhysicalName = CLBASE.template getMarkerNameByLocalBcId< CLT,2 >(); \
                    auto Expression = CLBASE.template getExprByLocalBcId<CLT,2>();\
                    EXPR;                                               \
                    break;}                                             \
                case 3 :  {                                             \
                    auto PhysicalName = CLBASE.template getMarkerNameByLocalBcId< CLT,3 >(); \
                    auto Expression = CLBASE.template getExprByLocalBcId<CLT,3>(); \
                    EXPR;                                               \
                    break;}                                             \
                case 4 :  {                                             \
                    auto PhysicalName = CLBASE.template getMarkerNameByLocalBcId< CLT,4 >(); \
                    auto Expression = CLBASE.template getExprByLocalBcId<CLT,4>();\
                    EXPR;                                               \
                    break;}                                             \
                }                                                       \
        }                                                               \
    /**/





    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    class cldataScal
    {
    public :
        static const bool is_scalar = true;
        static const bool is_vectorial = false;
        static const int nType = -1;
        static const int nName = -1;
        typedef vf::Expr< vf::Cst<double> > expression_type;

        std::string
        getMarkerName() const { return std::string("default-data-scalar"); }

        expression_type
        getExpr() const
        {
            return vf::cst(0.);
        }

    };

    //--------------------------------------------------------------------------------//

    class cldataVec
    {
    public :

        static const bool is_scalar = false;
        static const bool is_vectorial = true;
        static const int nType = -2;
        static const int nName = -2;
        typedef vf::Expr< vf::One<-1> > expression_type;

        std::string
        getMarkerName() const { return std::string("default-data-vectorial"); }

        expression_type
        getExpr() const
        {
            return vf::one();
        }

    };

    //--------------------------------------------------------------------------------//

    struct cldataGenBase
    {
        static const int scalaire = -1;
        static const int vectoriel = -2;

        template <int ShapeType>
        struct BcDataDefaultGenType
        {
            //BOOST_MPL_ASSERT_MSG( (ShapeType  == -1) || (ShapeType == -2), INVALID_DEFAULT_TYPE, ( mpl::int_<ShapeType> ) );
            typedef typename mpl::if_< boost::is_same<mpl::int_<ShapeType>,mpl::int_<scalaire> >,
                                       cldataScal,
                                       typename mpl::if_< boost::is_same<mpl::int_<ShapeType>,mpl::int_<vectoriel> >,
                                                          cldataVec,cldataVec >::type>::type type;
            //typedef typename type::expression_type expression_type;
        };

    };

    template <int TypeId>
    struct BcDataSpecific
    {
        static const int shape = cldataGenBase::template BcDataDefaultGenType<TypeId>::type::nType;
    };
    template <> struct BcDataSpecific<cl::dirichlet> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::dirichlet_vec> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::dirichlet_x> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::dirichlet_y> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::dirichlet_z> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::paroi_mobile> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::neumann_scal> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::neumann_vec> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::slip> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::pressure> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::fbm> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::fbm_dirichlet> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::fbm_dirichlet_vec> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::fbm_moving_wall> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::fluid_outlet> { static const int shape = cldataGenBase::scalaire; };
    template <> struct BcDataSpecific<cl::robin_vec> { static const int shape = cldataGenBase::vectoriel; };
    template <> struct BcDataSpecific<cl::follower_pressure> { static const int shape = cldataGenBase::scalaire; };


    template <int Type>
    class cldataGen : public cldataGenBase
    {
    public:
        static const int nShapeType = BcDataSpecific<Type>::shape;

        typedef typename cldataGenBase::BcDataDefaultGenType<nShapeType>::type type;
        static const int nName = type::nName;
        static const int nType = type::nType;
        typedef typename type::expression_type expression_type;

    };

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////


    class cldataNull
    {
    public:

        static const int value = 0;

        static const int nType = -1;
        static const int nName = -1;

        typedef vf::Expr< vf::Cst<double> > expression_type;

        expression_type
        getExpr() const { return vf::cst(0.); }

    };

    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    template<typename VectorBcDataType> class CLdataBase;

    template</*typename BcDataHeadType,*/typename... Args>
    class CLdataBase2 : public CLdataBase<boost::fusion::vector</*BcDataHeadType,*/Args...> >
    {
        typedef CLdataBase<boost::fusion::vector</*BcDataHeadType,*/Args...> > super_type;
    public :
        template<typename ...BcDataList>
        CLdataBase2( BcDataList... bcdatalist )
            :
            super_type( fusion::make_vector(bcdatalist...) )
            {}
        CLdataBase2(CLdataBase2 const & M) = default;
    };


    template<typename TTTT0>
    class CLdataBase// : public boost::fusion::vector<Args...>
    {
    public :
        typedef CLdataBase<TTTT0> self_type;
        //typedef boost::fusion::vector<TTTT0,Args...> mysuper_type;
        typedef TTTT0 mysuper_type;

        template<typename BCDataType>
        struct BcDataIsInvalid : mpl::less< mpl::int_< BCDataType::nType>, mpl::int_<0> >::type
        {};
        typedef typename mpl::remove_if< mysuper_type, BcDataIsInvalid<_> >::type mysuper_clean_type;

        template<int N>
        struct Basis
        {
            typedef typename mpl::if_< mpl::greater_equal< mpl::int_<N>, mpl::int_<0> >,
                                       typename mpl::if_< mpl::greater< mpl::int_<mpl::size<mysuper_type>::value>, mpl::int_<N> >,
                                                          typename mpl::at_c<mysuper_type,N>::type,
                                                          typename mpl::identity<cldataNull>::type >::type,
                                       typename cldataGenBase::BcDataDefaultGenType<N>::type >::type type;
            static const int nType = type::nType;
            static const int nName = type::nName;
        };
        template<typename BCDataType>
        struct ChangeTypeToSharedPtrType
        {
            typedef boost::shared_ptr<BCDataType> type;
        };
        typedef typename mpl::transform<mysuper_type, ChangeTypeToSharedPtrType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type container_type;

        template<typename BCDataType>
        struct ChangeBCDataTypeToIndices
        {
            typedef mpl::int_<BCDataType::nType> type;
        };
        typedef typename mpl::transform<mysuper_type, ChangeBCDataTypeToIndices<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type cldata_indices_type2;
        typedef typename mpl::remove_if< cldata_indices_type2, boost::is_same< _,mpl::int_<(-1)> > >::type cldata_indices_type;
        //static const bool check_cldata_indices = boost::is_same<cldata_indices_type,cldata_indices_type3>::value;
        //BOOST_MPL_ASSERT_MSG( check_cldata_indices,TYPE_NOT_EQUAL, ( mpl::size<cldata_indices_type>,mpl::size<cldata_indices_type2>,mpl::size<cldata_indices_type3> ) );


        template< int nType >
        struct hasThisType
        {
            static const bool value = mpl::contains<cldata_indices_type,mpl::int_<nType> >::type::value;
        };

        template <int Type,int numLoc>
        struct LocalToGlobalIndice
        {
            template <int nbF,int leN,int lecpt>
            struct ing2
            {
                static const int nbFind = nbF;
                static const int leNum = leN;
                static const int cpt = lecpt;
                static const int value = nbFind;
            };

            template<typename C,typename B/*,typename numLocBis*/>
            struct ingNext2
            {
                typedef typename mpl::if_< B ,
                                           ing2<C::nbFind+1,
                                                mpl::if_< mpl::equal_to< mpl::int_<C::nbFind>,mpl::int_<numLoc> >,
                                                          mpl::int_<C::cpt>,
                                                          mpl::int_<C::leNum> >::type::value,
                                               C::cpt+1 >,
                                           ing2<C::nbFind,C::leNum,C::cpt+1 > >::type type;
            };

            static const int value = mpl::fold< cldata_indices_type,
                                                ing2<0,cldataGen<Type>::nType/*0*/,0>,
                                                mpl::if_< mpl::equal_to<  mpl::_2 ,mpl::int_<Type> >,
                                                          ingNext2<mpl::_1,mpl::bool_<true> /*,mpl::int_<numLoc>*/ >,
                                                          ingNext2<mpl::_1,mpl::bool_<false> /*,mpl::int_<numLoc>*/ > >
                                                >::type::leNum;
        };

        template <int TheBcType,int TheName>
        struct KeyToGlobalIndice
        {
            template<typename BcDataType>
            struct SearchKey
            {
                typedef boost::is_same< fusion::vector< mpl::int_<TheBcType>, mpl::int_<TheName> >,
                                        fusion::vector< mpl::int_<BcDataType::nType>, mpl::int_<BcDataType::nName> > > type;
            };
            typedef typename mpl::find_if< mysuper_type, SearchKey<mpl::_1>  >::type find_key_type;
            static const int value = mpl::if_< boost::is_same<find_key_type,typename mpl::end<mysuper_type>::type >,
                                               mpl::int_<cldataGen<TheBcType>::nType>,
                                               mpl::int_<find_key_type::iter::pos::value> >::type::value;
        };

        ////////////////////////////////////////////////////////////////////////////////////////
        //constructors

        CLdataBase( container_type const& thecontainer )
            :
            M_datas( thecontainer)
            {}

        CLdataBase(CLdataBase const & M) = default;

        ////////////////////////////////////////////////////////////////////////////////////////


        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
        // hasThisType access
        static bool hasDirichlet() { return hasThisType<cl::dirichlet>::value; }
        static bool hasDirichletVec() { return hasThisType<cl::dirichlet_vec>::value; }
        static bool hasDirichletX() { return hasThisType<cl::dirichlet_x>::value; }
        static bool hasDirichletY() { return hasThisType<cl::dirichlet_y>::value; }
        static bool hasDirichletZ() { return hasThisType<cl::dirichlet_z>::value; }
        static bool hasMovingBoundary() { return hasThisType<cl::paroi_mobile>::value; }
        static bool hasNeumannScal() { return hasThisType<cl::neumann_scal>::value; }
        static bool hasNeumannVec() { return hasThisType<cl::neumann_vec>::value; }
        static bool hasSlipBoundary() { return hasThisType<cl::slip>::value; }
        static bool hasPressure() { return hasThisType<cl::pressure>::value; }
        static bool hasFbm() { return hasThisType<cl::fbm>::value; }
        static bool hasFbmDirichlet() { return hasThisType<cl::fbm_dirichlet>::value; }
        static bool hasFbmDirichletVec() { return hasThisType<cl::fbm_dirichlet_vec>::value; }
        static bool hasFbmMovingWall() { return hasThisType<cl::fbm_moving_wall>::value; }
        static bool hasFluidOutlet() { return hasThisType<cl::fluid_outlet>::value; }
        static bool hasRobinVec() { return hasThisType<cl::robin_vec>::value; }
        static bool hasFollowerPressure() { return hasThisType<cl::follower_pressure>::value; }

        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
        // objects access

        template<int N>
        boost::shared_ptr<typename Basis<N>::type>
        getData() const
        {
            return this->getData<N>( typename mpl::greater_equal< mpl::int_<Basis<N>::nType>, mpl::int_<0> >::type() );
        }
        template<int N>
        boost::shared_ptr<typename Basis<N>::type>
        getData(mpl::false_) const
        {
            //CHECK(false) << "not allowed";
            typedef typename Basis<N>::type type;
            return boost::shared_ptr<type>(new type());
        }
        template<int N>
        boost::shared_ptr<typename Basis<N>::type>
        getData(mpl::true_) const
        {
            return fusion::at_c<N>( M_datas );
        }
        //------------------------------------------------------------------------------------------------------------//
        template <int GlobalId>
        typename Basis< GlobalId >::type::expression_type
        getExpr() const
        {
            return this->getData<GlobalId>()->getExpr();
        }
        template <int BcTypeId,int MarkerId>
        typename Basis< KeyToGlobalIndice<BcTypeId,MarkerId>::value >::type::expression_type
        getExprByKey() const
        {
            return this->getExpr< KeyToGlobalIndice<BcTypeId,MarkerId>::value >();
        }
        template <int BcTypeId,int LocalId>
        typename Basis< LocalToGlobalIndice<BcTypeId,LocalId>::value >::type::expression_type
        getExprByLocalBcId() const
        {
            return this->getExpr< LocalToGlobalIndice<BcTypeId,LocalId>::value >();
        }
        //------------------------------------------------------------------------------------------------------------//
        template <int GlobalId>
        std::string
        getMarkerName() const
        {
            return this->getData<GlobalId>()->getMarkerName();
        }
        template <int BcTypeId,int MarkerId>
        std::string
        getMarkerNameByKey() const
        {
            return this->getMarkerName< KeyToGlobalIndice<BcTypeId,MarkerId>::value >();
        }
        template <int BcTypeId,int LocalId>
        std::string
        getMarkerNameByLocalBcId() const
        {
            return this->getMarkerName< LocalToGlobalIndice<BcTypeId,LocalId>::value >();
        }
        //------------------------------------------------------------------------------------------------------------//


        struct ForEachBCImplMarkerName
        {
            ForEachBCImplMarkerName( self_type const& bcDataBase )
                :
                M_bcDataBase( bcDataBase ),
                M_markerList()
                {}
            template<typename GlobalIndiceType>
            void operator()( GlobalIndiceType const& /*_globalIndices*/ ) const
            {
                std::string mark = M_bcDataBase.template getMarkerName<GlobalIndiceType::value>();
                M_markerList.push_back( mark );
            }
            std::list<std::string> markerList() const { return M_markerList; }
        private :
            self_type const& M_bcDataBase;
            mutable std::list<std::string> M_markerList;
        };

        template <int Type>
        std::list<std::string>
        getMarkerNameList() const
        {
            auto vectorIndiceByBc = this->globalIndicesByBc<Type>();
            ForEachBCImplMarkerName forEachFactory(*this);
            fusion::for_each( vectorIndiceByBc, forEachFactory );
            return forEachFactory.markerList();
        }
        //------------------------------------------------------------------------------------------------------------//
        template< int Type,int numLoc>
        int
        localToGlobal() const
        {
            return LocalToGlobalIndice<Type,numLoc>::value;
        }
        //------------------------------------------------------------------------------------------------------------//
        template <int T>
        int
        nbCL() const
        {
            return mpl::count_if< cldata_indices_type,
                                  boost::is_same< _ , mpl::int_<T> > >::type::value;
        }

        template<int TypeId>
        struct ForEachBcType
        {
            //typedef mpl::range_c<int,0,mpl::size<cldata_indices_type>::value > myrange_base_type;
            typedef mpl::range_c<int,0,mpl::size<mysuper_type>::value > myrange_base_type;

            template<typename RangeType>
            struct ContainerRangeType
            {
                typedef mpl::int_< RangeType::value > type;
            };
            typedef typename mpl::transform<myrange_base_type, ContainerRangeType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type myrange_type;

            template<typename RangeType>
            struct FilteredRangeType
            {
                typedef mpl::bool_<Basis<RangeType::value>::nType != TypeId> type;
            };
            typedef typename mpl::transform<myrange_type, FilteredRangeType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type myrange_filterd_type;

            template<typename RangeType>
            struct CleanRangeType
            {
                typedef mpl::bool_< mpl::at_c< myrange_filterd_type,RangeType::value >::type::value > type;
            };
            typedef typename mpl::remove_if< myrange_type, CleanRangeType<_> >::type myrange_clean_type;
            typedef myrange_clean_type type;
        };

        template<int TypeId>
        typename ForEachBcType<TypeId>::type
        globalIndicesByBc() const
        {
            return typename ForEachBcType<TypeId>::type();
        }
        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
        // operator+
        template<typename DataBaseType,int Nop=0>
        struct op_plus_return_type
        {
            typedef mpl::range_c<int,0,mpl::size<typename DataBaseType::mysuper_type>::value > myrange_base_type;

            template<typename RangeType>
            struct ContainerRangeType
            {
                typedef mpl::int_< RangeType::value > type;
            };
            typedef typename mpl::transform<myrange_base_type, ContainerRangeType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type myrange_type;

            template<typename RangeType>
            struct FilteredRangeType
            {
                typedef mpl::greater< mpl::int_<Nop>,RangeType > type;
            };
            typedef typename mpl::transform<myrange_type, FilteredRangeType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type myrange_filterd_type;

            template<typename RangeType>
            struct CleanRangeType
            {
                typedef mpl::bool_< mpl::at_c< myrange_filterd_type,RangeType::value >::type::value > type;
            };
            typedef typename mpl::remove_if< myrange_type, CleanRangeType<_> >::type myrange_clean_type;

            template<typename RangeType>
            struct CreateBcDataType
            {
                typedef typename mpl::at_c<typename DataBaseType::mysuper_type,RangeType::value >::type type;
            };
            typedef typename mpl::transform<myrange_clean_type, CreateBcDataType<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type container_to_add_type;

            typedef typename mpl::copy< container_to_add_type,
                                        mpl::back_inserter< mysuper_type > >::type concatenated_container_type;

            typedef CLdataBase<concatenated_container_type> type2;

            typedef typename mpl::push_back< mysuper_type, typename DataBaseType::template Basis<Nop>::type >::type add_head_only_container_type;
            typedef CLdataBase<add_head_only_container_type> add_head_only_type;

        };




        template<typename BcDataTypeAddedType>
        typename op_plus_return_type<BcDataTypeAddedType>::type2
        operator+( BcDataTypeAddedType const & dat)
        {
            return this->opPlus<0>(dat);
        }

        template<int Nop,typename BcDataTypeAddedType>
        typename op_plus_return_type<BcDataTypeAddedType,Nop>::type2
        opPlus(BcDataTypeAddedType const & opdat)
        {
            return opPlus<Nop>(opdat, mpl::bool_< (mpl::size<typename BcDataTypeAddedType::mysuper_type>::value == Nop) >() );
        }

        template<int Nop, typename BcDataTypeAddedType>
        typename op_plus_return_type<BcDataTypeAddedType,Nop>::type2
        opPlus(BcDataTypeAddedType const & opdat, mpl::false_)
        {
            auto newseq = fusion::push_back(this->M_datas,opdat.template getData<0>() );
            typename op_plus_return_type<BcDataTypeAddedType,Nop>::add_head_only_type hola( newseq );
            return hola.template opPlus<Nop+1>( opdat );
        }
        template<int Nop,typename BcDataTypeAddedType>
        container_type //typename op_plus_return_type<BcDataTypeAddedType,Nop>::type2
        opPlus(BcDataTypeAddedType const & opdat, mpl::true_)
        {
            return M_datas;
        }
        // fin operator+
        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
    private :
        container_type M_datas;

    };



        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////

    template <int clType,int clName,typename ExprType>
    class CLdata  : public boost::enable_shared_from_this< CLdata<clType,clName,ExprType> >
    {
    public :

        typedef CLdata<clType,clName,ExprType> self_type;

        static const int nType = clType;
        static const int nName = clName;

        typedef ExprType expression_type;
        typedef boost::shared_ptr<expression_type> expression_ptrtype;

        CLdata()
        {
            //M_clString = boost::make_tuple(cl::cl_type_string[clType],cl::cl_marker_string[clName]);
            M_clString=cl::cl_marker_string[clName];
        }

        CLdata(ExprType const & expr)
            :
            M_clExpr(new expression_type(expr))
        {
            //M_clString=boost::make_tuple(cl::cl_type_string[clType],cl::cl_marker_string[clName]);
            M_clString=cl::cl_marker_string[clName];
        }

        CLdata(const CLdata & M)
            :
            M_clExpr(M.M_clExpr),
            M_clString(M.M_clString)
        {}

#if 0
        CLdataBase2<self_type>
        convertToBase()
        {
            return CLdataBase2<self_type>(this->shared_from_this());
        }

        std::string
        getTypeName() const { return boost::get<0>(M_clString); }
#endif
        std::string
        //getMarkerName() const { return boost::get<1>(M_clString); }
        getMarkerName() const { return M_clString; }

        expression_type const &
        getExpr() const { return *M_clExpr; }

    private :

        expression_ptrtype M_clExpr;

        //boost::tuple<std::string,std::string> M_clString;
        std::string M_clString;

    };








    template <int cltype,int clname,typename ExprType>
    CLdataBase2<CLdata<cltype,clname,ExprType> >
    addCL(ExprType expr)
    {

        typedef CLdata<cltype,clname,ExprType> cl_un_type;
        boost::shared_ptr<CLdata<cltype,clname,ExprType> > clbis(new CLdata<cltype,clname,ExprType>(expr));
        auto temp = boost::shared_ptr< CLdataBase2< cl_un_type> >(new CLdataBase2< cl_un_type>(clbis->shared_from_this()));
        return *temp;
    }

    template <int cltype,int clname>
    CLdataBase2<CLdata<cltype,clname,vf::Expr< vf::Cst<double> > > >
    addCL()
    {
        typedef vf::Expr< vf::Cst<double> > ExprType;

        typedef CLdata<cltype,clname,ExprType> cl_un_type;
        boost::shared_ptr<CLdata<cltype,clname,ExprType> > clbis(new CLdata<cltype,clname,ExprType>(vf::cst(0.)));
        //boost::shared_ptr<CLdata<cltype,clname,ExprType> > clbis(new CLdata<cltype,clname,ExprType>(vf::Px()));
        auto temp = boost::shared_ptr< CLdataBase2< cl_un_type> >(new CLdataBase2< cl_un_type>(clbis->shared_from_this()));
        return *temp;
    }



} // end namespace Feel


#endif /* __CLDATA_H */
