/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Stephane Veys <stephane.veys@imag.fr>
 Date: 2013-02-22

 Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)
 Copyright (C) 2011-present Feel++ Consortium
 
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
   \file model.hpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-02-22
*/
#ifndef ModelCrbBase_H
#define ModelCrbBase_H 1

#include <feel/feelcrb/crbmodeldb.hpp>
#include <feel/feelcrb/eim.hpp>
#include <feel/feelcrb/deim.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feelvf/vf.hpp>
//#include<boost/tokenizer.hpp>
#include<boost/regex.hpp>
namespace Feel
{

enum {
    /** TimeIndependent */
    TimeIndependent=0,
    /** TimeDependent */
    TimeDependent = 0x1,
    /**  */
    Linear = 0,
    /**  */
    NonLinear = 0x2,
    /** Coercive PDE */
    Coercive = 0,
    /** Inf-Sup PDE */
    InfSup = 0x4,
    /** Monolithic */
    Mono = 0,
    /** By block */
    UseBlock=0x8
};

class ModelCrbBaseBase : public CommObject
{
public:
    using super = CommObject;
    ModelCrbBaseBase() : super( Environment::worldCommPtr() ) {}
    explicit ModelCrbBaseBase( worldcomm_ptr_t const& w ) : super( w ) {}
    ModelCrbBaseBase( ModelCrbBaseBase const& ) = default;
    ModelCrbBaseBase( ModelCrbBaseBase && ) = default;
    ModelCrbBaseBase& operator=( ModelCrbBaseBase const& ) = default;
    ModelCrbBaseBase& operator=( ModelCrbBaseBase && ) = default;
    virtual ~ModelCrbBaseBase() = default;
};

class ParameterDefinitionBase
{
public :
    typedef ParameterSpace<1> parameterspace_type ;
};

class FunctionSpaceDefinitionBase
{
public :
    /*mesh*/
    typedef Simplex<1> entity_type ;
    typedef Mesh<entity_type > mesh_type ;

    /*basis*/
    typedef bases< Lagrange<1,Scalar> > basis_type ;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;

    static const bool is_time_dependent=false;
    static const bool is_linear=true;
};

template <typename ParameterDefinition, typename FunctionSpaceDefinition>
class EimDefinitionBase
{

public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename mpl::if_<is_shared_ptr<FunctionSpaceDefinition>,
                              mpl::identity<typename FunctionSpaceDefinition::element_type>,
                              mpl::identity<FunctionSpaceDefinition>>::type::type::space_type space_type;
    //typedef typename FunctionSpaceDefinition::space_type space_type;

    /* EIM */
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fun_type ;
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fund_type ;

};



template <typename ParameterDefinition=ParameterDefinitionBase,
          typename FunctionSpaceDefinition=FunctionSpaceDefinitionBase,
          int _Options = 0,
          typename EimDefinition = EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition>
          >
class ModelCrbBase :
        public ModelCrbBaseBase,
        public std::enable_shared_from_this< ModelCrbBase<ParameterDefinition,FunctionSpaceDefinition,_Options,EimDefinition> >
{

public :
    using super = ModelCrbBaseBase;
    typedef ModelCrbBase self_type;
    typedef typename EimDefinition::fun_type fun_type;
    typedef typename EimDefinition::fund_type fund_type;

    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef typename parameterspace_type::element_type parameter_type;

    typedef typename mpl::if_<is_shared_ptr<FunctionSpaceDefinition>,
                              mpl::identity<typename FunctionSpaceDefinition::element_type>,
                              mpl::identity<FunctionSpaceDefinition>>::type::type::space_type space_type;
    typedef space_type functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef functionspace_ptrtype space_ptrtype;

    typedef typename functionspace_type::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename functionspace_type::basis_type basis_type;
    typedef typename functionspace_type::value_type value_type;

    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    /*reduced basis space*/
    typedef ReducedBasisSpace<space_type> rbfunctionspace_type;
    typedef std::shared_ptr< rbfunctionspace_type > rbfunctionspace_ptrtype;
    typedef typename rbfunctionspace_type::ctxrbset_type rbfunctionspace_context_type;
    typedef typename rbfunctionspace_type::ctxrbset_ptrtype rbfunctionspace_context_ptrtype;

    /*backend*/
    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;


    typedef std::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    typedef std::shared_ptr<fund_type> fund_ptrtype;
    typedef std::vector<fund_ptrtype> funsd_type;

    typedef Eigen::VectorXd vectorN_type;

    typedef OperatorLinear< space_type , space_type > operator_type;
    typedef std::shared_ptr<operator_type> operator_ptrtype;

    typedef OperatorLinearComposite< space_type , space_type > operatorcomposite_type;
    typedef std::shared_ptr<operatorcomposite_type> operatorcomposite_ptrtype;

    typedef FsFunctionalLinearComposite< space_type > functionalcomposite_type;
    typedef std::shared_ptr<functionalcomposite_type> functionalcomposite_ptrtype;

    typedef FsFunctionalLinear< space_type > functional_type;
    typedef std::shared_ptr<functional_type> functional_ptrtype;


    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef DEIMBase<parameterspace_type,space_type,vector_type> deim_type;
    typedef std::shared_ptr<deim_type> deim_ptrtype;
    typedef DEIMBase<parameterspace_type,space_type,sparse_matrix_type> mdeim_type;
    typedef std::shared_ptr<mdeim_type> mdeim_ptrtype;

    typedef std::vector<deim_ptrtype> deim_vector_type;
    typedef std::vector<mdeim_ptrtype> mdeim_vector_type;

    //static const uint16_type ParameterSpaceDimension = ParameterDefinition::ParameterSpaceDimension ;
    typedef std::vector< std::vector< double > > beta_vector_type;
    typedef std::vector< double > beta_vector_light_type;

    typedef vf::detail::BilinearForm<functionspace_type, functionspace_type,VectorUblas<value_type>> form2_type;
    typedef vf::detail::LinearForm<functionspace_type,vector_type,vector_type> form1_type;

    typedef Pchv_type<mesh_type,1> displacement_space_type;
    typedef std::shared_ptr<displacement_space_type> displacement_space_ptrtype;
    typedef typename displacement_space_type::element_type displacement_field_type;
    typedef typename std::shared_ptr<displacement_field_type> displacement_field_ptrtype;

#if 0
    static const bool is_time_dependent = FunctionSpaceDefinition::is_time_dependent;
    static const bool is_linear = FunctionSpaceDefinition::is_linear;
#else
    static const bool is_time_dependent = ((_Options&TimeDependent)==TimeDependent);
    //static const bool is_linear = ((_Options&Linear)==Linear);
    static const bool is_linear = !((_Options&NonLinear)==NonLinear);
    static const bool by_block = (_Options&UseBlock)==UseBlock;
#endif
    static const int Options = _Options;





    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Mq
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Aq
                                   std::vector< std::vector<std::vector<vector_ptrtype> > >//Fq
                                   >,
                               boost::tuple<
                                   std::vector< std::vector<sparse_matrix_ptrtype> >,//Aq
                                   std::vector< std::vector<std::vector<vector_ptrtype> > >//Fq
                                   >
                               >::type affine_decomposition_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   std::vector< sparse_matrix_ptrtype >,//Mq
                                   std::vector< sparse_matrix_ptrtype >,//Aq
                                   std::vector< std::vector< vector_ptrtype > >//Fq
                                   >,
                               boost::tuple<
                                   std::vector< sparse_matrix_ptrtype >,//Aq
                                   std::vector< std::vector< vector_ptrtype > >//Fq
                                   >
                               >::type affine_decomposition_light_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   sparse_matrix_ptrtype, sparse_matrix_ptrtype, std::vector<vector_ptrtype>
                                   >,
                               boost::tuple<
                                   sparse_matrix_ptrtype, std::vector<vector_ptrtype>
                                   >
                               >::type monolithic_type;


    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   std::map<int,double>, std::map<int,double>, std::vector< std::map<int,double> >
                                   >,
                               boost::tuple<
                                   std::map<int,double>, std::vector< std::map<int,double> >
                                   >
                               >::type eim_interpolation_error_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >,
                               boost::tuple<
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   >
                               >::type betaqm_type;

    typedef typename mpl::if_< mpl::bool_< is_time_dependent >,
                               boost::tuple<
                                   beta_vector_light_type,
                                   beta_vector_light_type,
                                   std::vector<beta_vector_light_type>
                                   >,
                               boost::tuple<
                                   beta_vector_light_type,
                                   std::vector<beta_vector_light_type>
                                   >
                               >::type betaq_type;

    typedef std::vector< boost::tuple< sparse_matrix_ptrtype , std::string > > lhs_light_type;
    typedef std::shared_ptr<lhs_light_type> lhs_light_ptrtype;
    typedef std::vector< boost::tuple<  vector_ptrtype , std::string > > rhs_light_type;
    typedef std::shared_ptr<rhs_light_type> rhs_light_ptrtype ;
    typedef std::vector< boost::tuple<  std::vector< vector_ptrtype > , std::string > > outputs_light_type;
    typedef std::shared_ptr<outputs_light_type> outputs_light_ptrtype ;

    typedef Bdf<space_type>  bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;

    ModelCrbBase() = delete;
    ModelCrbBase( std::string const& name,  worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(), std::string const& prefix = "" )
        :
        ModelCrbBase( name, Environment::randomUUID( true ), worldComm, prefix )
        {}
    ModelCrbBase( std::string const& name, uuids::uuid const& uid, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(), std::string const& prefix = "" )
        :
        super( worldComm ),
        M_prefix(prefix),
        Dmu( parameterspace_type::New( 0,worldComm ) ),
        XN( new rbfunctionspace_type( worldComm ) ),
        M_backend( backend() ),
        M_crbModelDb( name, uid ),
        M_is_initialized( false ),
        M_has_displacement_field( false ),
        M_rebuildDb(boption(_prefix=M_prefix,_name="crb.rebuild-database")),
        M_dbLoad(ioption(_prefix=M_prefix, _name="crb.db.load" )),
        M_dbFilename(soption(_prefix=M_prefix, _name="crb.db.filename")),
        M_dbId(soption(_prefix=M_prefix,_name="crb.db.id")),
        M_dbUpdate(ioption(_prefix=M_prefix, _name="crb.db.update" )),
        M_serEimFreq(ioption(_prefix=M_prefix,_name = "ser.eim-frequency")),
        M_serRbFreq(ioption(_prefix=M_prefix,_name = "ser.rb-frequency")),
        M_serErrorEstimation(boption(_prefix=M_prefix,_name="ser.error-estimation")),
        M_crbUseNewton(boption(_prefix=M_prefix,_name="crb.use-newton")),
        M_isOnlineModel(false)
    {

        if ( !M_rebuildDb )
        {
            switch ( M_dbLoad )
            {
            case 0 :
                M_crbModelDb.updateIdFromDBFilename( M_dbFilename );
                break;
            case 1:
                M_crbModelDb.updateIdFromDBLast( crb::last::created );
                break;
            case 2:
                M_crbModelDb.updateIdFromDBLast( crb::last::modified );
                break;
            case 3:
                M_crbModelDb.updateIdFromId( M_dbId );
                break;
            }
        }
        else
        {
            switch ( M_dbUpdate )
            {
            case 0 :
                M_crbModelDb.updateIdFromDBFilename( M_dbFilename );
                break;
            case 1:
                M_crbModelDb.updateIdFromDBLast( crb::last::created );
                break;
            case 2:
                M_crbModelDb.updateIdFromDBLast( crb::last::modified );
                break;
            case 3:
                M_crbModelDb.updateIdFromId( M_dbId );
                break;
            default:
                // don't do anything and let the system pick up a new unique id
                break;
            }
        }

        if ( this->worldComm().isMasterRank() )
        {
            fs::path modeldir = M_crbModelDb.dbRepository();
            Feel::cout << "Model repository: " << modeldir.string() << std::endl;
            if ( !fs::exists( modeldir ) )
                fs::create_directories( modeldir );
            boost::property_tree::ptree ptree;
            ptree.add( "name", this->modelName() );
            ptree.add( "uuid", uuids::to_string(this->uuid()) );
            std::string jsonpath = (modeldir / M_crbModelDb.jsonFilename()).string();
            write_json( jsonpath, ptree );
        }
        this->worldComm().barrier();

    }

    virtual ~ModelCrbBase() {}

    /**
     * \return the name of the model
     */
    std::string const& modelName() const { return M_crbModelDb.name(); }

    /**
     * set the model name
     */
    void setModelName( std::string const& name ) { M_crbModelDb.setName( name ); }

    std::string const& prefix() const { return M_prefix; }

    void setPrefix( std::string const& prefix ) { M_prefix = prefix; }

    /**
     * Define the model as an online (sequential) model
     * used for the computation of coefficient during the online phase
     */
    void setModelOnlineDeim( std::string const& name )
    {
        M_backend = backend( _name=name, _worldcomm=this->subWorldCommSeqPtr() );
    }

    //! functions call by deim to init specific part of the model when online.
    virtual void initOnlineModel()
    {}

    void setOnlineModel() { M_isOnlineModel = true; }
    bool isOnlineModel() const { return M_isOnlineModel; }

    //!
    //! unique id for CRB Model
    //!
    uuids::uuid  uuid() const { return M_crbModelDb.uuid(); }

    //!
    //! set uuid for CRB Model
    //! @warning be extra careful here, \c setId should be called before any
    //! CRB type object is created because they use the id
    //!
    void setId( uuids::uuid const& i ) { M_crbModelDb.setId( i ); }

    /**
     * return directory path where symbolic expression (ginac) are built
     * empty means current path
     */
    std::string const& symbolicExpressionBuildDir() const { return M_symbolicExpressionBuildDir; }
    /**
     * set directory path where symbolic expression (ginac) are built
     */
    void setSymbolicExpressionBuildDir( std::string const& dir ) { M_symbolicExpressionBuildDir=dir; }

    /**
     * list of files to copy in database
     */
    std::map<std::string,std::string> const& additionalModelFiles() const { return M_additionalModelFiles; }

    /**
     * add file to copy in database
     */
    void addModelFile( std::string const& key, std::string const& filename )
        {
            M_additionalModelFiles[key] = filename;
        }

    /**
     * return true if model has registered a file with this key
     */
    bool hasModelFile( std::string const& key ) const
        {
            return ( M_additionalModelFiles.find( key ) != M_additionalModelFiles.end() );
        }


    void addLhs( boost::tuple< form2_type, std::string > const & tuple )
    {
        M_Aq.push_back( tuple.template get<0>().matrixPtr() );
        M_betaAqString.push_back( tuple.template get<1>() );
    }
    /*
     * return the left hand side terms from the affine decomposition
     * that is to say bilinear forms Aq and beta coefficients associated
     */
    void addLhs( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple )
    {
        M_Aq.push_back( tuple.template get<0>() );
        M_betaAqString.push_back( tuple.template get<1>() );
    }

    /*
     * return the terms from the affine decomposition linked to the mass matrix
     * that is to say bilinear forms Mq and beta coefficients associated
     */
    void addMass( boost::tuple< form2_type , std::string > const & tuple )
    {
        M_Mq.push_back( tuple.template get<0>().matrixPtr() );
        M_betaMqString.push_back( tuple.template get<1>() );
    }

    /*
     * return the terms from the affine decomposition linked to the mass matrix
     * that is to say bilinear forms Mq and beta coefficients associated
     */
    void addMass( boost::tuple< sparse_matrix_ptrtype , std::string > const & tuple )
    {
        M_Mq.push_back( tuple.template get<0>() );
        M_betaMqString.push_back( tuple.template get<1>() );
    }

    /*
     * return the right hand side terms from the affine decomposition
     * that is to say linear forms Fq[0] and beta coefficients associated
     */
    void addRhs( boost::tuple<  form1_type , std::string > const & tuple )
    {
        if( M_Fq.size() == 0 )
        {
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[0].push_back( tuple.template get<0>().vectorPtr() );
        M_betaFqString[0].push_back( tuple.template get<1>() );
    }

    /*
     * return the right hand side terms from the affine decomposition
     * that is to say linear forms Fq[0] and beta coefficients associated
     */
    void addRhs( boost::tuple<  vector_ptrtype , std::string > const & tuple )
    {
        if( M_Fq.size() == 0 )
        {
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[0].push_back( tuple.template get<0>() );
        M_betaFqString[0].push_back( tuple.template get<1>() );
    }

    /*
     * return the terms from the affine decomposition linked to output
     * that is to say linear forms and beta coefficients associated
     */
    void addOutput( boost::tuple<  form1_type , std::string > const & tuple )
    {
        int size = M_Fq.size();
        if( size == 0 )
        {
            //there is no rhs or output yet
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[1].push_back( tuple.template get<0>().vectorPtr() );
        M_betaFqString[1].push_back( tuple.template get<1>() );
    }
    /*
     * return the terms from the affine decomposition linked to output
     * that is to say linear forms and beta coefficients associated
     */
    void addOutput( boost::tuple<  vector_ptrtype , std::string > const & tuple )
    {
        int size = M_Fq.size();
        if( size == 0 )
        {
            //there is no rhs or output yet
            M_Fq.resize(2);
            M_betaFqString.resize(2);
        }
        M_Fq[1].push_back( tuple.template get<0>() );
        M_betaFqString[1].push_back( tuple.template get<1>() );
    }

    void addEnergyMatrix( form2_type const & f )
    {
        M_energy_matrix = f.matrixPtr() ;
    }
    void addEnergyMatrix( sparse_matrix_ptrtype const & matrix )
    {
        M_energy_matrix = matrix ;
    }

    void addMassMatrix( form2_type const & f )
    {
        M_mass_matrix = f.matrixPtr() ;
    }
    void addMassMatrix( sparse_matrix_ptrtype const & matrix )
    {
        M_mass_matrix = matrix ;
    }

    void setInitialized( const bool & b )
    {
        M_is_initialized = b ;
    }

    bool isInitialized()
    {
        return M_is_initialized;
    }

    std::map<std::string,std::vector<boost::any> >
    rbSpaceContextEim() const
    {
        std::map<std::string,std::vector<boost::any> > res;
        for ( int k=0;k<M_funs.size();++k )
        {
            //boost::any const& ctxRbAny = M_funs[k]->rbSpaceContext();
            if ( M_funs[k]->hasRbSpaceContext() )
                res[M_funs[k]->name()].push_back( M_funs[k]->rbSpaceContext() );
            if ( M_funs[k]->hasRbSpaceContext2() )
                res[M_funs[k]->name()].push_back( M_funs[k]->rbSpaceContext2() );
        }
        for ( int k=0;k<M_funs_d.size();++k )
        {
            if ( M_funs_d[k]->hasRbSpaceContext() )
                res[M_funs_d[k]->name()].push_back( M_funs_d[k]->rbSpaceContext() );
            if ( M_funs_d[k]->hasRbSpaceContext2() )
                res[M_funs_d[k]->name()].push_back( M_funs_d[k]->rbSpaceContext2() );
        }
        return res;
    }

#if 0
    std::map<std::string,rbfunctionspace_context_ptrtype>
    rbSpaceContextEim() const
    {
        std::map<std::string,rbfunctionspace_context_ptrtype> res;
        for ( int k=0;k<M_funs.size();++k )
        {
            boost::any const& ctxRbAny = M_funs[k]->rbSpaceContext();
            if ( !boost::any_cast<rbfunctionspace_context_ptrtype>( &ctxRbAny ) )
                continue;
            auto ctxRb = boost::any_cast<rbfunctionspace_context_ptrtype>( ctxRbAny );
            res[M_funs[k]->name()] = ctxRb;
        }
        for ( int k=0;k<M_funs_d.size();++k )
        {
            boost::any const& ctxRbAny = M_funs_d[k]->rbSpaceContext();
            if ( !boost::any_cast<rbfunctionspace_context_ptrtype>( &ctxRbAny ) )
                continue;
            auto ctxRb = boost::any_cast<rbfunctionspace_context_ptrtype>( ctxRbAny );
            res[M_funs_d[k]->name()] = ctxRb;
        }
        return res;
    }
#endif
    void
    setRbSpaceContextEim( std::map<std::string,std::vector<boost::any> > const& rbCtx )
    {
        for ( int k=0;k<M_funs.size();++k )
        {
            std::string const& eimName = M_funs[k]->name();
            auto itFindRbCtx = rbCtx.find( eimName );
            if ( itFindRbCtx == rbCtx.end() )
                continue;
            auto const& rbCtxVec = itFindRbCtx->second;
            if ( rbCtxVec.size() > 0 )
                M_funs[k]->setRbSpaceContext( rbCtxVec[0] );
            if ( rbCtxVec.size() > 1 )
                M_funs[k]->setRbSpaceContext2( rbCtxVec[1] );
        }
        for ( int k=0;k<M_funs_d.size();++k )
        {
            std::string const& eimName = M_funs_d[k]->name();
            auto itFindRbCtx = rbCtx.find( eimName );
            if ( itFindRbCtx == rbCtx.end() )
                continue;
            auto const& rbCtxVec = itFindRbCtx->second;
            if ( rbCtxVec.size() > 0 )
                M_funs_d[k]->setRbSpaceContext( rbCtxVec[0] );
            if ( rbCtxVec.size() > 1 )
                M_funs_d[k]->setRbSpaceContext2( rbCtxVec[1] );
        }
    }
#if 0
    void
    setRbSpaceContextEim( std::map<std::string,rbfunctionspace_context_ptrtype> const& rbCtx )
    {
        for ( int k=0;k<M_funs.size();++k )
        {
            std::string const& eimName = M_funs[k]->name();
            auto itFindRbCtx = rbCtx.find( eimName );
            if ( itFindRbCtx == rbCtx.end() )
                continue;
            rbfunctionspace_context_ptrtype rbCtx2 = itFindRbCtx->second;
            rbCtx2->setRbFunctionSpace( this->rBFunctionSpace() );
            M_funs[k]->setRbSpaceContext( rbCtx2 );
        }
        for ( int k=0;k<M_funs_d.size();++k )
        {
            std::string const& eimName = M_funs_d[k]->name();
            auto itFindRbCtx = rbCtx.find( eimName );
            if ( itFindRbCtx == rbCtx.end() )
                continue;
            rbfunctionspace_context_ptrtype rbCtx2 = itFindRbCtx->second;
            rbCtx2->setRbFunctionSpace( this->rBFunctionSpace() );
            M_funs_d[k]->setRbSpaceContext( rbCtx2 );
        }
    }
#endif
    void updateRbSpaceContextEim()
    {
        for ( int k=0;k<M_funs.size();++k )
            M_funs[k]->updateRbSpaceContext( this->rBFunctionSpace() );
        for ( int k=0;k<M_funs_d.size();++k )
            M_funs_d[k]->updateRbSpaceContext( this->rBFunctionSpace() );
    }

    void addEim( fun_ptrtype const& fun )
    {
        M_funs.push_back( fun );
    }

    void addEimDiscontinuous( fund_ptrtype const& fund)
    {
        M_funs_d.push_back( fund );
    }

    virtual funs_type scalarContinuousEim () const
    {
        return M_funs;
    }

    virtual funsd_type scalarDiscontinuousEim () const
    {
        return M_funs_d;
    }

    void clearScalarContinuousEim ()
    {
        M_funs.clear();
    }

    void clearScalarDiscontinuousEim ()
    {
        M_funs_d.clear();
    }

    void addDeim( deim_ptrtype const& deim )
    {
        M_deims.push_back( deim );
    }
    void addMdeim( mdeim_ptrtype const& mdeim )
    {
        M_mdeims.push_back( mdeim);
    }

    deim_ptrtype deim( int const& i=0 ) const
    {
        return M_deims[i];
    }

    mdeim_ptrtype mdeim( int const& i=0 ) const
    {
        return M_mdeims[i];
    }

    bool hasDeim() const
    {
        return M_deims.size() || M_mdeims.size();
    }

    deim_vector_type deimVector() const
    {
        return M_deims;
    }
    mdeim_vector_type mdeimVector() const
    {
        return M_mdeims;
    }

    virtual vector_ptrtype assembleForDEIM( parameter_type const& mu, int const& tag )
    {
        return nullptr;
    }

    virtual vector_ptrtype assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
    {
        return nullptr;
    }


    virtual sparse_matrix_ptrtype assembleForMDEIM( parameter_type const& mu, int const& tag )
    {
        return nullptr;
    }

    virtual sparse_matrix_ptrtype assembleForMDEIMnl( parameter_type const& mu, element_type const& u, int const& tag )
    {
        return nullptr;
    }


    /**
     * \brief update model description in property_tree
     * \param ptree to update
     */
    void updatePropertyTree( boost::property_tree::ptree & ptree ) const
    {
        ptree.add( "name", this->modelName());

        boost::property_tree::ptree ptreeParameterSpace;
        this->parameterSpace()->updatePropertyTree( ptreeParameterSpace );
        ptree.add_child( "parameter_space", ptreeParameterSpace );

        boost::property_tree::ptree ptreeModelFiles;
        for ( auto const& filenamePair : this->additionalModelFiles() )
        {
            boost::property_tree::ptree ptreeModelFile;
            std::string const& key = filenamePair.first;
            std::string const& filename = filenamePair.second;
            ptreeModelFile.add("filename",filename );
            ptreeModelFiles.add_child( key, ptreeModelFile );
        }
        ptree.add_child( "additional-model-files", ptreeModelFiles );

        boost::property_tree::ptree ptreeEim;
        for ( auto const& eimObject : this->scalarContinuousEim() )
        {
            boost::property_tree::ptree ptreeEimObject;
            // ptreeEimObject.add("database-filename", (fs::path(eimObject->dbDirectory())/fs::path(eimObject->dbFilename())).string() );
            ptreeEimObject.add("database-filename", eimObject->dbRelativePath() );
            ptreeEim.add_child( eimObject->name(), ptreeEimObject );
        }
        for ( auto const& eimObject : this->scalarDiscontinuousEim() )
        {
            boost::property_tree::ptree ptreeEimObject;
            // ptreeEimObject.add("database-filename", (fs::path(eimObject->dbDirectory())/fs::path(eimObject->dbFilename())).string() );
            ptreeEimObject.add("database-filename", eimObject->dbRelativePath() );
            ptreeEim.add_child( eimObject->name(), ptreeEimObject );
        }
        ptree.add_child( "eim", ptreeEim );

        boost::property_tree::ptree ptreeDeim;
        for ( auto const& deimObject : this->deimVector() )
        {
            boost::property_tree::ptree ptreeDeimObject;
            ptreeDeimObject.add("database-filename", deimObject->dbRelativePath() );
            ptreeDeim.add_child( deimObject->name(), ptreeDeimObject );
        }
        for ( auto const& mdeimObject : this->mdeimVector() )
        {
            boost::property_tree::ptree ptreeMdeimObject;
            ptreeMdeimObject.add("database-filename", mdeimObject->dbRelativePath() );
            ptreeDeim.add_child( mdeimObject->name(), ptreeMdeimObject );
        }
        ptree.add_child( "deim", ptreeDeim );


        boost::property_tree::ptree ptree_betaFq;
        for (int k=0;k<M_betaFq.size();++k )
        {
            boost::property_tree::ptree ptree_betaFq_value;
            ptree_betaFq_value.put( "",M_betaFq[k].size() );
            ptree_betaFq.push_back( std::make_pair("", ptree_betaFq_value) );
        }
        boost::property_tree::ptree ptree_betaAqm;
        for (int k=0;k<M_betaAqm.size();++k )
        {
            boost::property_tree::ptree ptree_betaAqm_value;
            ptree_betaAqm_value.put( "",M_betaAqm[k].size() );
            ptree_betaAqm.push_back( std::make_pair("", ptree_betaAqm_value) );
        }
        boost::property_tree::ptree ptree_betaMqm;
        for (int k=0;k<M_betaMqm.size();++k )
        {
            boost::property_tree::ptree ptree_betaMqm_value;
            ptree_betaMqm_value.put( "",M_betaMqm[k].size() );
            ptree_betaMqm.push_back( std::make_pair("", ptree_betaMqm_value) );
        }
        boost::property_tree::ptree ptree_betaFqm;
        for (int k=0;k<M_betaFqm.size();++k )
        {
            boost::property_tree::ptree ptree_betaFqm_sub;
            for (int k2=0;k2<M_betaFqm[k].size();++k2 )
            {
                boost::property_tree::ptree ptree_betaFqm_sub_value;
                ptree_betaFqm_sub_value.put( "",M_betaFqm[k][k2].size() );
                ptree_betaFqm_sub.push_back( std::make_pair("", ptree_betaFqm_sub_value) );
            }
            ptree_betaFqm.push_back( std::make_pair("", ptree_betaFqm_sub) );
        }

        boost::property_tree::ptree ptreeAffineDecomposition;
        ptreeAffineDecomposition.add( "betaAq", M_betaAq.size() );
        ptreeAffineDecomposition.add( "betaMq", M_betaMq.size() );
        ptreeAffineDecomposition.add_child( "betaFq", ptree_betaFq );
        ptreeAffineDecomposition.add_child( "betaAqm", ptree_betaAqm );
        ptreeAffineDecomposition.add_child( "betaMqm", ptree_betaMqm );
        ptreeAffineDecomposition.add_child( "betaFqm", ptree_betaFqm );
        ptree.add_child( "affine-decomposition", ptreeAffineDecomposition );

        boost::property_tree::ptree ptreeSpecificityOfModel;
        this->updateSpecificityModel( ptreeSpecificityOfModel );
        ptree.add_child( "specifity-of-model", ptreeSpecificityOfModel );
    }
    /**
     * \brief load CrbModel from json
     * \param input json filename
     */
    void loadJson( std::string const& filename, std::string const& childname = "" )
    {
        if ( !fs::exists( filename ) )
        {
            LOG(INFO) << "Could not find " << filename << std::endl;
            return;
        }

        std::string dbDir = fs::path( filename ).parent_path().string();

        auto json_str_wo_comments = removeComments(readFromFile(filename));
        //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

        boost::property_tree::ptree ptree;
        std::istringstream istr( json_str_wo_comments );
        boost::property_tree::read_json( istr, ptree );
        if ( childname.empty() )
            this->setup( ptree, dbDir );
        else
        {
            auto const& ptreeChild = ptree.get_child( childname );
            this->setup( ptreeChild, dbDir );
        }
    }
    void setup( boost::property_tree::ptree const& ptree, std::string const& dbDir )
    {
        auto const& ptreeParameterSpace = ptree.get_child( "parameter_space" );
        Dmu->setup( ptreeParameterSpace );

        auto ptreeModelFiles = ptree.get_child_optional( "additional-model-files" );
        if ( ptreeModelFiles )
            for ( auto const& ptreeModelFilePair : *ptreeModelFiles/*ptree.get_child( "additional-model-files" )*/ )
            {
                std::string const& key = ptreeModelFilePair.first;
                auto const& ptreeModelFile = ptreeModelFilePair.second;
                std::string filenameAdded  = ptreeModelFile.template get<std::string>( "filename" );
                fs::path filenameAddedPath = fs::path( filenameAdded );
                if ( !dbDir.empty() && filenameAddedPath.is_relative() )
                    filenameAddedPath = fs::path(dbDir) / filenameAddedPath;

                // filenameAdded = (fs::path(dbDir)/fs::path(filenameAdded).filename()).string();

                // std::cout << "additional-model-files : key=" << key << " filename=" << filenameAddedPath.string() << "\n";
                this->addModelFile( key, filenameAddedPath.string()/*filenameAdded*/ );
            }

        auto ptreeAffineDecomposition = ptree.get_child_optional( "affine-decomposition" );
        if ( ptreeAffineDecomposition )
        {
            int sizeBetaAq = ptreeAffineDecomposition->template get<int>( "betaAq" );
            M_betaAq.resize( sizeBetaAq );
            int sizeBetaMq = ptreeAffineDecomposition->template get<int>( "betaMq" );
            M_betaMq.resize( sizeBetaMq );
            int sizeBetaFq = ptreeAffineDecomposition->get_child("betaFq").size();
            M_betaFq.resize( sizeBetaFq );
            int k=0;
            for( auto const& item : ptreeAffineDecomposition->get_child("betaFq") )
            {
                int sizeSubBetaFq = item.second.template get_value<int>();
                M_betaFq[k].resize( sizeSubBetaFq );
                ++k;
            }
            int sizeBetaAqm = ptreeAffineDecomposition->get_child("betaAqm").size();
            M_betaAqm.resize( sizeBetaAqm );
            k=0;
            for( auto const& item : ptreeAffineDecomposition->get_child("betaAqm") )
            {
                int sizeSubBetaAqm = item.second.template get_value<int>();
                M_betaAqm[k].resize( sizeSubBetaAqm );
                ++k;
            }
            int sizeBetaMqm = ptreeAffineDecomposition->get_child("betaMqm").size();
            M_betaMqm.resize( sizeBetaMqm );
            k=0;
            for( auto const& item : ptreeAffineDecomposition->get_child("betaMqm") )
            {
                int sizeSubBetaMqm = item.second.template get_value<int>();
                M_betaMqm[k].resize( sizeSubBetaMqm );
                ++k;
            }
            int sizeBetaFqm = ptreeAffineDecomposition->get_child("betaFqm").size();
            M_betaFqm.resize( sizeBetaFqm );
            k=0;
            for( auto const& item : ptreeAffineDecomposition->get_child("betaFqm") )
            {
                int sizeSubBetaFqm = item.second.get_child("").size();
                M_betaFqm[k].resize( sizeSubBetaFqm );
                int k2=0;
                for ( auto const& item2 : item.second.get_child("") )
                {
                    int thesize = item2.second.template get_value<int>();
                    M_betaFqm[k][k2].resize( thesize );
                    ++k2;
                }
                ++k;
            }

        }


        this->setupSpecificityModel( ptree, dbDir );

        XN->setFunctionSpace( this->functionSpace() );
        //XN->setModel( this->shared_from_this() );

        this->setInitialized( true );
    }

    virtual void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) {}
    virtual void updateSpecificityModel( boost::property_tree::ptree & ptree ) const {}

    virtual void initModel() = 0;

    virtual void assemble()
    {
        if( (M_serEimFreq != 0) || (M_serRbFreq != 0) )
        {
            Feel::cout << "************************************************************************ \n"
                       << "** SER method is asked without implementation of assemble() function. ** \n"
                       << "** assemble() function is needed to update affine decomposition during \n"
                       << "** the simultaneous basis functions build ** \n"
                       << "************************************************************************" << std::endl;
        }
    };

    virtual eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN )
    {
        return eimInterpolationErrorEstimation( mu, uN,  mpl::bool_<is_time_dependent>() );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN , mpl::bool_<true> )
    {
        return boost::make_tuple( M_eim_error_mq , M_eim_error_aq, M_eim_error_fq);
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN , mpl::bool_<false> )
    {
        return boost::make_tuple( M_eim_error_aq, M_eim_error_fq);
    }

    virtual eim_interpolation_error_type eimInterpolationErrorEstimation( )
    {
        return eimInterpolationErrorEstimation( mpl::bool_<is_time_dependent>() );
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( mpl::bool_<true> )
    {
        return boost::make_tuple( M_eim_error_mq , M_eim_error_aq, M_eim_error_fq);
    }
    eim_interpolation_error_type eimInterpolationErrorEstimation( mpl::bool_<false> )
    {
        return boost::make_tuple( M_eim_error_aq, M_eim_error_fq);
    }


    virtual std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA()
    {
        if( M_Aqm.size() == 0 && this->worldComm().isMasterRank() )
        {
            std::cout<<"************************************************************************"<<std::endl;
            std::cout<<"** It seems that you are using operators free and you don't have      **"<<std::endl;
            std::cout<<"** implemented computeLinearDecompositionA() to have a linear         **"<<std::endl;
            std::cout<<"** decomposition of the bilinear form.                                **"<<std::endl;
            std::cout<<"** It will be used to compute norm of the error during CRB convergence**"<<std::endl;
            std::cout<<"************************************************************************"<<std::endl;
        }
        if( M_linearAqm.size() == 0 )
        {
            CHECK( false )<<"The linear part of a() for models using EIM must be filled ! You have to implement computeLinearDecompositionA() \n";
        }
        return M_linearAqm;
    }

    /*
     * the user has to provide the affine decomposition
     * he has two choices : use vectors/matrices or operator free
     * if he uses operators free he must implement :
     * - operatorCompositeM (only if the model is time dependent)
     * - operatorCompositeA
     * - functionalCompositeF
     * else, if he uses vectors/matrices (classical way) he must implement
     * - computeAffineDecomposition
     */

    virtual operatorcomposite_ptrtype operatorCompositeA()
    {
        return M_compositeA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeM()
    {
        return M_compositeM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeF()
    {
        return M_compositeF;
    }

    virtual operatorcomposite_ptrtype operatorCompositeLightA()
    {
        return M_compositeA;
    }
    virtual operatorcomposite_ptrtype operatorCompositeLightM()
    {
        return M_compositeM;
    }
    virtual std::vector< functionalcomposite_ptrtype > functionalCompositeLightF()
    {
        return M_compositeF;
    }

    virtual affine_decomposition_type computeAffineDecomposition()
    {
        return computeAffineDecomposition( mpl::bool_< is_time_dependent >() );
    }

    affine_decomposition_type computeAffineDecomposition( mpl::bool_<true> )
    {
        return boost::make_tuple( M_Mqm , M_Aqm , M_Fqm );
    }

    affine_decomposition_type computeAffineDecomposition( mpl::bool_<false> )
    {
        return boost::make_tuple( M_Aqm , M_Fqm );
    }

    virtual affine_decomposition_type computePicardAffineDecomposition()
    {
        if( this->worldComm().isMasterRank() && M_serErrorEstimation && M_crbUseNewton )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** Use of SER error estimation with newton needs     **"<<std::endl;
            std::cout<<"** computePicardAffineDecomposition(...) **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        return computeAffineDecomposition();
    }

    virtual affine_decomposition_light_type computeAffineDecompositionLight()
    {
        return computeAffineDecompositionLight( mpl::bool_< is_time_dependent >() );
    }
    affine_decomposition_light_type computeAffineDecompositionLight( mpl::bool_<true> )
    {
        return boost::make_tuple( M_Mq , M_Aq , M_Fq );
    }
    affine_decomposition_light_type computeAffineDecompositionLight( mpl::bool_<false> )
    {
        return boost::make_tuple( M_Aq , M_Fq );
    }

    virtual monolithic_type computeMonolithicFormulation( parameter_type const& mu )
    {
        return computeMonolithicFormulation( mu , mpl::bool_< is_time_dependent >() );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu, mpl::bool_<true> )
    {
        return boost::make_tuple( M_monoM , M_monoA , M_monoF );
    }
    monolithic_type computeMonolithicFormulation( parameter_type const& mu, mpl::bool_<false> )
    {
        return boost::make_tuple( M_monoA , M_monoF );
    }

    virtual monolithic_type computeMonolithicFormulationU( parameter_type const& mu, element_type const& u )
    {
        return computeMonolithicFormulationU( mu , u, mpl::bool_< is_time_dependent >() );
    }
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu, element_type const& u, mpl::bool_<true> )
    {
        return boost::make_tuple( M_monoM , M_monoA , M_monoF );
    }
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu, element_type const& u, mpl::bool_<false> )
    {
        return boost::make_tuple( M_monoA , M_monoF );
    }

    virtual beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu ,  double time=0 )
    {
        return computeBetaLinearDecompositionA( mu, mpl::bool_< is_time_dependent >(), time );
    }
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, mpl::bool_<true>, double time )
    {
        auto tuple = computeBetaQm( mu , time );
        return tuple.template get<1>();
    }
    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu, mpl::bool_<false>, double time )
    {
        auto tuple = computeBetaQm( mu , time );
        return tuple.template get<0>();
    }

    // Default updateResidual / updateJacobian functions
    virtual bool updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm){ return false; }
    virtual bool updateJacobian(element_type const& X, std::vector< std::vector<sparse_matrix_ptrtype> >& Jqm){ return false; }

    virtual betaq_type computeBetaQ( parameter_type const& mu ,  double time , bool only_terms_time_dependent=false )
    {
        return computeBetaQ( mu, mpl::bool_< is_time_dependent >(), time );
    }
    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<true>, double time , bool only_terms_time_dependent=false )
    {
        int sizeA=M_ginacAq.size();
        if( sizeA > 0 )
        {
            //if the user doesn't implement computeBetaQ function
            //and so gives affine decomposition terms using addLhs ect...
            //then we consider we have 2 outputs
            int nboutputs=2;
            int sizeM=M_ginacMq.size();
            int sizeF0=M_ginacFq[0].size();
            int sizeF1=M_ginacFq[1].size();
            int musize=mu.size();
            std::string symbol;

            std::map<std::string,double> map_symbols;
            for(int i=0; i<musize; i++)
            {
                map_symbols.insert( std::pair<std::string, double>(Dmu->parameterName(i), mu(i)));
            }

            M_betaAq.resize( sizeA );
            for(int q=0; q<sizeA; q++)
            {
                double coeff = M_ginacAq[q].evaluate(map_symbols)(0,0);
                M_betaAq[q]=coeff;
            }
            M_betaMq.resize( sizeM );
            for(int q=0; q<sizeM; q++)
            {
                double coeff = M_ginacMq[q].evaluate( map_symbols )(0,0);
                M_betaMq[q]=coeff;
            }
            M_betaFq.resize( nboutputs );
            for(int output=0; output<nboutputs;output++)
            {
                int size=M_betaFqString[output].size();
                M_betaFq[output].resize( size );
                for(int q=0; q<size; q++)
                {
                    double coeff = M_ginacFq[output][q].evaluate( map_symbols )(0,0);
                    M_betaFq[output][q]=coeff;
                }
            }
        }

        return boost::make_tuple( M_betaMq, M_betaAq, M_betaFq );
    }
    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<false>, double time , bool only_terms_time_dependent=false)
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** Error ! You want to access to computeBetaQ ( mu , time) but   **"<<std::endl;
            std::cout<<"** your model is not time-dependent !                            **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        bool go=false;
        CHECK( go );
        return boost::make_tuple( M_betaAq, M_betaFq );
    }

    //for steady models, only mu is needed to compute beta coefficients
    virtual betaqm_type computeBetaQm( parameter_type const& mu )
    {
        return computeBetaQm( mu, mpl::bool_< is_time_dependent >() );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<true>  )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** Error ! You want to access to computeBetaQm( mu ) wherease    **"<<std::endl;
            std::cout<<"** your model is time-dependent !                                **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        bool go=false;
        CHECK( go );
        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false> )
    {
        return boost::make_tuple( M_betaAqm, M_betaFqm );
    }

    virtual betaqm_type computePicardBetaQm( parameter_type const& mu )
    {
        if( this->worldComm().isMasterRank() && M_serErrorEstimation && M_crbUseNewton )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** Use of SER error estimation with newton needs"<<std::endl;
            std::cout<<"** computePicardBetaQm(...) function **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        return computeBetaQm( mu );
    }

    virtual betaq_type computeBetaQ( parameter_type const& mu )
    {
        return computeBetaQ( mu, mpl::bool_< is_time_dependent >() );
    }

    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<true>  )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** Error ! You want to access to computeBetaQ ( mu ) wherease    **"<<std::endl;
            std::cout<<"** your model is time-dependent !                                **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        bool go=false;
        CHECK( go );
        return boost::make_tuple( M_betaMq, M_betaAq, M_betaFq );
    }
    betaq_type computeBetaQ( parameter_type const& mu , mpl::bool_<false> )
    {
        int sizeA=M_ginacAq.size();

        if( sizeA > 0 )
        {
            //if the user doesn't implement computeBetaQ function
            //and so gives affne decomposition terms using addLhs ect...
            //then we consider we have 2 outputs
            int nboutputs=2;
            int sizeM=M_ginacMq.size();
            int sizeF0=M_ginacFq[0].size();
            int sizeF1=M_ginacFq[1].size();
            int musize=mu.size();
            std::string symbol;

            std::map<std::string,double> map_symbols;
            for(int i=0; i<musize; i++)
            {
                map_symbols.insert( std::pair<std::string, double>(Dmu->parameterName(i), mu(i)));
            }

            M_betaAq.resize( sizeA );
            for(int q=0; q<sizeA; q++)
            {
                double coeff = M_ginacAq[q].evaluate(map_symbols)(0,0);
                M_betaAq[q]=coeff;
            }
            M_betaFq.resize( nboutputs );
            for(int output=0; output<nboutputs;output++)
            {
                int size=M_betaFqString[output].size();
                M_betaFq[output].resize( size );
                for(int q=0; q<size; q++)
                {
                    double coeff = M_ginacFq[output][q].evaluate( map_symbols )(0,0);
                    M_betaFq[output][q]=coeff;
                }
            }
        }
        return boost::make_tuple( M_betaAq, M_betaFq );
    }


    virtual betaqm_type computeBetaQm( parameter_type const& mu ,  double time , bool only_terms_time_dependent=false)
    {
        return computeBetaQm( mu, mpl::bool_< is_time_dependent >(), time , only_terms_time_dependent );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<true>, double time , bool only_terms_time_dependent=false )
    {
        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }
    betaqm_type computeBetaQm( parameter_type const& mu , mpl::bool_<false>, double time , bool only_terms_time_dependent=false )
    {
        return boost::make_tuple( M_betaAqm, M_betaFqm );
    }

    virtual betaqm_type computePicardBetaQm( parameter_type const& mu ,  double time , bool only_terms_time_dependent=false)
    {
        if( this->worldComm().isMasterRank() && M_serErrorEstimation && M_crbUseNewton )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** Use of SER error estimation with newton needs"<<std::endl;
            std::cout<<"** computePicardBetaQm(...) function **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        return computeBetaQm( mu, time , only_terms_time_dependent );
    }



    virtual std::vector<double> computeBetaTri( parameter_type const& mu )
        {
            return M_betaTri;
        }

    virtual sparse_matrix_ptrtype assembleTrilinearJ( element_type const& u, int const& q )
        {
            CHECK(false) << "You have to implement assembleTrilinearJ to use trilinear crb models\n";
            sparse_matrix_ptrtype M;
            return M;
        }
    virtual vector_ptrtype assembleTrilinearR( element_type const& u, int const& q )
        {
            CHECK(false) << "You have to implement assembleTrilinearR to use trilinear crb models\n";
            vector_ptrtype R;
            return R;
        }
    virtual sparse_matrix_ptrtype assembleTrilinearOperator( element_type const& u, int const& q )
        {
            CHECK(false) << "You have to implement assembleTrilinearOperator to use trilinear crb models\n";
            sparse_matrix_ptrtype M;
            return M;
        }

    virtual int sizeOfBilinearJ() const
        { return -1; }

    virtual int sizeOfLinearR() const
        { return -1; }

    virtual bool hasZeroMeanPressure()
        { return false; }

    void buildGinacBetaExpressions( parameter_type const& mu )
    {
        //not that the parameter mu is here to indicates
        //the dimension of the parameter space
        //and to construct symbols
        int musize = mu.size();
        for(int i=0; i<musize; i++)
        {
            M_symbols_vec.push_back( Dmu->parameterName(i) );
        }
        fs::path dir( Environment::expand( M_symbolicExpressionBuildDir ) );
        if( M_betaAqString.size() > 0 )
        {
            int size = M_betaAqString.size();
            for(int q=0; q<size; q++)
            {
                std::string filename = ( boost::format("GinacA%1%") %q ).str();
                if ( !M_symbolicExpressionBuildDir.empty() )
                    filename = ( dir / fs::path( filename ) ).string();
                M_ginacAq.push_back( expr( M_betaAqString[q], Symbols( M_symbols_vec ) , filename ) );
            }
        }
        if( M_betaMqString.size() > 0 )
        {
            int size = M_betaMqString.size();
            for(int q=0; q<size; q++)
            {
                std::string filename = ( boost::format("GinacM%1%") %q ).str();
                if ( !M_symbolicExpressionBuildDir.empty() )
                    filename = ( dir / fs::path( filename ) ).string();
                M_ginacMq.push_back( expr( M_betaMqString[q], Symbols( M_symbols_vec ) , filename ) );
            }
        }
        if( M_betaFqString.size() > 0 )
        {
            int nboutputs = M_betaFqString.size();
            M_ginacFq.resize( nboutputs );
            for(int output=0; output<nboutputs; output++)
            {
                int size = M_betaFqString[output].size();
                for(int q=0; q<size; q++)
                {
                    std::string filename = ( boost::format("GinacF%1%.%2%") %output %q ).str();
                    if ( !M_symbolicExpressionBuildDir.empty() )
                        filename = ( dir / fs::path( filename ) ).string();
                    M_ginacFq[output].push_back( expr( M_betaFqString[output][q], Symbols( M_symbols_vec ) , filename ) );
                }
            }
        }
    }


    //this function is not called bdf() to not interfere with bdf constructor
    virtual bdf_ptrtype bdfModel()
    {
        if( is_time_dependent )
        {
            if( this->worldComm().isMasterRank() )
            {
                std::cout<<"*******************************************************************"<<std::endl;
                std::cout<<"** Error ! You have implemented a transient problem but you      **"<<std::endl;
                std::cout<<"** forgot to implement bdfModel() function that returns your bdf **"<<std::endl;
                std::cout<<"*******************************************************************"<<std::endl;
            }
            bool go=false;
            CHECK( go );
        }
        bdf_ptrtype dummy_bdf;
        return dummy_bdf;
    }

    //initialize the field for transient problems
    virtual void initializationField( element_ptrtype& initial_field,parameter_type const& mu )
    {
        initial_field->setZero();
    }


    //for linear models, beta coefficients don't depend on solution u
    //so the user doesn't have to specify this function
    virtual betaqm_type computeBetaQm( element_type const& u, parameter_type const& mu ,  double time , bool only_time_dependent_terms=false )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** You are using the function computeBetaQm( u , mu , time ) but **"<<std::endl;
            std::cout<<"** your model has only implemented computeBetaQm( mu , time )    **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        betaqm_type dummy_beta_coeff;
        return dummy_beta_coeff;
    }
    virtual betaqm_type computeBetaQm( element_type const& u , parameter_type const& mu )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** You are using the function computeBetaQm( u , mu ) but     **"<<std::endl;
            std::cout<<"** your model has only implemented computeBetaQm( mu )        **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        betaqm_type dummy_beta_coeff;
        return dummy_beta_coeff;
    }
    virtual betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu , double time , bool only_time_dependent_terms=false )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** You are using the function computeBetaQm( u , mu , time ) but **"<<std::endl;
            std::cout<<"** your model doesnt implement this function                     **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        betaqm_type dummy_beta_coeff;
        return dummy_beta_coeff;
    }
    virtual betaqm_type computeBetaQm( vectorN_type const& urb, parameter_type const& mu )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"*******************************************************************"<<std::endl;
            std::cout<<"** You are using the function computeBetaQm( urb , mu ) but      **"<<std::endl;
            std::cout<<"** your model doesnt implement this function                     **"<<std::endl;
            std::cout<<"*******************************************************************"<<std::endl;
        }
        betaqm_type dummy_beta_coeff;
        return dummy_beta_coeff;
    }

    virtual betaqm_type computePicardBetaQm( element_type const& u , parameter_type const& mu )
    {
        if( this->worldComm().isMasterRank() && M_serErrorEstimation && M_crbUseNewton )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** Use of SER error estimation with newton needs"<<std::endl;
            std::cout<<"** computePicardBetaQm(...) function **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        return this->computeBetaQm( u, mu );
    }
    virtual betaqm_type computePicardBetaQm( element_type const& u, parameter_type const& mu ,  double time , bool only_time_dependent_terms=false )
    {
        if( this->worldComm().isMasterRank() && M_serErrorEstimation && M_crbUseNewton )
        {
            std::cout<<"****************************************************************"<<std::endl;
            std::cout<<"** Use of SER error estimation with newton needs"<<std::endl;
            std::cout<<"** computePicardBetaQm(...) function **"<<std::endl;
            std::cout<<"****************************************************************"<<std::endl;
        }
        return computeBetaQm( u, mu, time, only_time_dependent_terms );
    }



    /**
     * By default reference parameters are min parameters
     * If the user needs to specify them
     * then this function should returns true
     */
    virtual bool referenceParametersGivenByUser()
    {
        return false;
    }
    virtual parameter_type refParameter()
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"**************************************************************************************************"<<std::endl;
            std::cout<<"** You want to specify reference parameters because referenceParametersGivenByUser returns true **"<<std::endl;
            std::cout<<"** your must impelement the function refParameter() !                                           **"<<std::endl;
            std::cout<<"**************************************************************************************************"<<std::endl;
        }
        bool go=false;
        CHECK( go );
        parameter_type muref;
        return muref;
    }

    //for linear steady models, mass matrix does not exist
    //non-linear steady models need mass matrix for the initial guess
    //this class can provide the function operatorCompositeM() to ensure compilation
    //but we need to know if the model can provide an operator composite for mass matrix
    //because non-linear problems have to do these operators
    //because CRBModel is not able to construct such operators in the general case.
    //Elements of function space are members of CRBModel
    //then for composite spaces, we need a view of these elements
    //BUT we can't have a view of an element of a non-composite space
    //this function returns true if the model provides an operator composite for mass matrix
    //false by default
    virtual bool constructOperatorCompositeM()
    {
        return false;
    }



    /**
     * note about the initial guess :
     * by default, if the model doesn't give an initial guess
     * then the initial guess is zero
     */
    virtual int QInitialGuess() const
    {
        return 1;
    }

    virtual int mMaxInitialGuess( int q ) const
    {
        if( q == 0 )
            return 1;
        else
            throw std::logic_error( "[ModelCrbBase::mMaxInitialGuess(q)] ERROR wrong index q, should be 0" );
    }

    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu )
    {
        beta_vector_type beta;
        beta.resize(1); //q=1
        beta[0].resize(1); //m=1
        beta[0][0]=0;
        return beta;
    }



    /**
     * solve the model for a given parameter \p mu
     * linear models don't need to provide this fuction
     * but non-linear : yes
     * and also those using CRBTrilinear
     */
    virtual std::pair<element_type,bool> safeSolve( parameter_type const& mu )
        {
            element_type u = solve( mu );
            return std::make_pair( u, true );
        }
    virtual element_type solve( parameter_type const& mu )
    {
        element_type solution;
        return solution;
    }


    /**
     * inner product
     */
    virtual sparse_matrix_ptrtype energyMatrix ()
    {
        double norm = M_energy_matrix->l1Norm();
        CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
        return M_energy_matrix;
    }

    virtual sparse_matrix_ptrtype energyMatrix () const
    {
        double norm = M_energy_matrix->l1Norm();
        CHECK( norm > 0 )<<"The energy matrix has not be filled !\n";
        return M_energy_matrix;
    }

    /**
     * inner product for mass matrix
     * Transient models need to implement these functions.
     */
    virtual sparse_matrix_ptrtype const& massMatrix () const
    {
        double norm = M_mass_matrix->l1Norm();
        CHECK( norm > 0 )<<"The mass matrix has not be filled !\n";
        return M_mass_matrix;
    }
    virtual sparse_matrix_ptrtype massMatrix ()
    {
        double norm = M_mass_matrix->l1Norm();
        CHECK( norm > 0 )<<"The mass matrix has not be filled !\n";
        return M_mass_matrix;
    }
#if 0
    virtual sparse_matrix_ptrtype const& innerProductForPod () const
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForPod function");
        return M;
    }
    virtual sparse_matrix_ptrtype innerProductForPod ()
    {
        throw std::logic_error("Your model is time-dependant so you MUST implement innerProductForPod function");
        return M;
    }
#endif

    /*
     * If the model is nonlinear and then need to implement an initial guess
     * it has to implement it !
     * Here is just a code to compile if model is a linear model.
     */
    virtual std::vector< std::vector< element_ptrtype > >
    computeInitialGuessAffineDecomposition()
    {
        throw std::logic_error("Your model is non linear so you MUST implement computeInitialGuessAffineDecomposition function");
        std::vector< std::vector<element_ptrtype> > q;
        return q;
    }


    void partitionMesh( std::string mshfile , std::string target , int dimension, int order )
    {
        int N = this->worldComm().globalSize();

        if( this->worldComm().isMasterRank() )
            std::cout<<"[ModelCrbBase] generate target file : "<<target<<" from "<<mshfile<<std::endl;

            Gmsh gmsh( dimension,
                       order,
                       Environment::worldCommPtr() );
           gmsh.setNumberOfPartitions( N );
           gmsh.rebuildPartitionMsh( mshfile /*mesh with 1 partition*/, target /*mesh with N partitions*/ );
    }


    /**
     * compute statistics on vectors
     * arguments : vectors of double and associated names
     * results : statistics on data
     */
    vectorN_type computeStatistics( Eigen::VectorXd const& vector , std::string const& name )
    {
        double min=0,max=0,mean=0,mean1=0,mean2=0,standard_deviation=0,variance=0;
        Eigen::MatrixXf::Index index;
        Eigen::VectorXd square;

        std::ofstream file;
        std::string filename="OnlineStatistics-"+name;

        file.open( filename , std::ios::app );

        if( vector.size() > 0 )
        {
            bool force = boption("eim.use-dimension-max-functions");
            int Neim=0;
            if( force )
                Neim = ioption("eim.dimension-max");

            int N = vector.size();

            min = vector.minCoeff(&index);
            max = vector.maxCoeff(&index);
            mean = vector.mean();
            mean1 = mean * mean;
            square  = vector.array().pow(2);
            mean2 = square.mean();
            standard_deviation = math::sqrt( mean2 - mean1 );

            if( this->worldComm().isMasterRank() )
            {
                if( force )
                {
                    std::cout<<"statistics  for "<<name<<" (  was called "<< N << " times with "<<Neim<<" basis functions )"<<std::endl;
                    file <<" statistics  for "<<name<<" (  was called "<< N << " times with "<<Neim<<" basis functions )\n";
                }
                else
                {
                    std::cout<<" statistics  for "<<name<<" (  was called "<< N << " times )"<<std::endl;
                    file <<" statistics  for "<<name<<" (  was called "<< N << " times )\n";
                }
                std::cout<<"min : "<<min<<" - max : "<<max<<" mean : "<<mean<<" standard deviation : "<<standard_deviation<<"  (see "<<filename<<")"<<std::endl;
                file<<"min : "<<min<<" - max : "<<max<<" mean : "<<mean<<" standard deviation : "<<standard_deviation<<"\n";
            }
        }
        vectorN_type result(4);
        result(0)=min;
        result(1)=max;
        result(2)=mean;
        result(3)=standard_deviation;
        return result;
    }


    void writeConvergenceStatistics( std::vector< vectorN_type > const& vector, std::string const& filename , std::string const& extra="")
    {
        if( this->worldComm().isMasterRank() )
        {
            Eigen::MatrixXf::Index index_max;
            Eigen::MatrixXf::Index index_min;

            double totaltime=0;
            double globaltotaltime=0;
            std::ofstream file;
            file.open( filename );
            if( extra == "totaltime" )
                file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\t"<< "Total time" << "\n";
            else
                file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\t" << "Mean" << "\t" << "Variance" << "\n";
            int Nmax = vector.size();
            std::vector<double> nbruns( Nmax );
            for(int n=0; n<Nmax; n++)
            {
                totaltime=0;
                double variance=0;
                int sampling_size = vector[n].size();
                double mean=vector[n].mean();
                double max = vector[n].maxCoeff(&index_max);
                double min = vector[n].minCoeff(&index_min);
                for(int i=0; i<sampling_size; i++)
                {
                    variance += (vector[n](i) - mean)*(vector[n](i) - mean)/sampling_size;
                    totaltime+= vector[n](i);
                }
                globaltotaltime += totaltime;
                if( extra == "totaltime" )
                    file <<n+1<<"\t"<<min<<"\t"<<max<<"\t"<<mean<<"\t"<<variance<<"\t"<<totaltime<<"\n";
                else
                    file <<n+1<<"\t"<<min<<"\t"<<max<<"\t"<<mean<<"\t"<<variance<<"\n";
                nbruns[n]=vector[n].size();
            }

            if( extra == "totaltime" )
            {
                file << "#global total time : "<<globaltotaltime<<"\n";
            }

            //write information about number of runs
            for(int n=0; n<Nmax; n++)
            {
                file <<"#N = "<<n<<" -- number of runs : "<<nbruns[n]<<"\n";
            }

        }
    }

    /**
    * \param : vector1
    * \param : vector2
    * look for the index minIdx for which we have min( vector1 )
    * then compute vector2( minIdx) / vector1( minIdx )
    * and same thing with max : i.e. vector2( maxIdx ) / vector1( maxIdx )
    * usefull to compute error estimation efficiency associated to
    * min/max errors
    **/
    void writeVectorsExtremumsRatio( std::vector< vectorN_type > const& error, std::vector< vectorN_type > const& estimated, std::string const& filename )
    {
        if( this->worldComm().isMasterRank() )
        {

            Eigen::MatrixXf::Index index_max;
            Eigen::MatrixXf::Index index_min;

            std::ofstream file;
            file.open( filename );
            file << "NbBasis" << "\t" << "Min" << "\t" << "Max" << "\n";
            int Nmax = error.size();
            int estimatedsize = estimated.size();
            CHECK( Nmax == estimatedsize ) << "vectors have not the same size ! v1.size() : "<<error.size()<<" and v2.size() : "<<estimated.size()<<"\n";
            for(int n=0; n<Nmax; n++)
            {
                double min_error = error[n].maxCoeff(&index_min);
                double max_error = error[n].minCoeff(&index_max);
                double min_estimated = estimated[n](index_min);
                double max_estimated = estimated[n](index_max);
                double min_ratio = min_estimated / min_error;
                double max_ratio = max_estimated / max_error;
                file << n+1 <<"\t"<< min_ratio <<" \t" << max_ratio<< "\n";
            }
        }
    }

    void readConvergenceDataFromFile( std::vector< vectorN_type > & vector, std::string const& filename )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::vector< std::vector< double > > tmpvector;
            std::ifstream file ( filename );
            //std::stringstream file ( filename );
            std::string str;
            int N,i;
            int Nmax=0;
            double value;
            boost::regex re("[+-]?([[:digit:]]*\\.)?[[:digit:]]+([eE][+-]?[[:digit:]]+)?");
            if( file )
            {
                //first, determine max elements of the RB
                //because we could have performed several runs
                //and the size of the RB can vary between two runs
                std::string _s;
                while( std::getline(file, _s) )
                {
                    std::vector< double > _v;
                    boost::sregex_iterator it(_s.begin(),_s.end(),re); //parse the list of numbers

                    std::string first=it->str();
                    N=std::atoi(first.c_str());
                    it++;
                    if( N > Nmax )
                        Nmax = N;

                    boost::sregex_iterator j;
                    i=0;
                    _v.resize(N);
                    for(; it!=j; ++it)
                    {
                        if( (*it)[0].matched )
                        {
                            std::string s = it->str();
                            auto value = std::atof(s.c_str());
                            if(value != std::numeric_limits<double>::min())
                            {
                                _v[i]=value;
                            }
                        }
                        else
                            throw std::logic_error( "[ModelCrbBase::fillVectorFromFile] ERROR : Cannot read value" );
                        i++;
                    }
                    //CHECK(i == N);
                    tmpvector.push_back(_v);
                }
                vector.resize( Nmax );
            }
            else
            {
                std::cout<<"The file "<<filename<<" was not found "<<std::endl;
                throw std::logic_error( "[ModelCrbBase::fillVectorFromFile] ERROR loading the file " );
            }
            file.close();

            //now copy std::vector into eigen vector
            int nbvalues = tmpvector.size();
            for(int n=0; n<Nmax; n++)
            {
                vector[n].resize(nbvalues);
                for(int i=0; i<nbvalues; i++)
                {
                    vector[n](i) = tmpvector[i][n];
                }
            }

        }//master proc
    }//end of function


    /*
     * if the model has geometric parameters then the mesh can be adapted to
     * the current parameter to visualize solution field
     */
    void adaptMesh( parameter_type const& mu ){ /*by default nothing to be done*/ ; }

    /*
     * \param components_vary : vector of indices for components vary
     * \param extremums : vector containing min and max parameters valuers
     * \param cuttings : vector containing the cutting in each direction + time initial and time final and time step used
     * \param time_vary : (bool) the time vary if true
     */
    gmsh_ptrtype createStructuredGrid( std::vector<int> components_vary, std::vector<parameter_type> extremums, std::vector<int> cuttings,
                                       std::vector<double> time_cuttings, bool time_vary)
    {
        auto min=extremums[0];
        auto max=extremums[1];
        double min0 = min(components_vary[0]);
        double min1 = min(components_vary[1]);
        double max0 = max(components_vary[0]);
        double max1 = max(components_vary[1]);
        double Ti=time_cuttings[0];
        double Tf=time_cuttings[1];
        double dt=time_cuttings[2];
        int nb0=cuttings[0];
        int nb1=cuttings[1];
        if( time_vary )
        {
            nb0=(Tf-Ti)/dt;
            min0=Ti;
            max0=Tf;
        }
        gmsh_ptrtype gmshp( new Gmsh );
        std::ostringstream ostr;

        //we want that each cell created here will be a cell of the mesh
        //so we take a large hsize
        int p=1;
        ostr <<"min0 = "<<min0<<";\n"
             <<"max0 = "<<max0<<";\n"
             <<"min1 = "<<min1<<";\n"
             <<"max1 = "<<max1<<";\n"
             <<"Point (1) = { min0, min1, 1, hsize };\n"
             <<"Point (2) = { max0, min1, 1, hsize };\n"
             <<"Point (3) = { max0, max1, 1, hsize };\n"
             <<"Point (4) = { min0, max1, 1, hsize };\n"
             <<"Line (1) = { 1 , 2 };\n"
             <<"Line (2) = { 2 , 3 };\n"
             <<"Line (3) = { 3 , 4 };\n"
             <<"Line (4) = { 4 , 1 };\n"
             <<"Line Loop (1) = {1,2,3,4};\n"
             <<"Plane Surface (100) = {1};\n"
             <<"Transfinite Line{1,-3} = "<<nb0<<";\n"
             <<"Transfinite Line{2,-4} = "<<nb1<<";\n"
             <<"Transfinite Surface{100} = {1,2,3,4};\n"
             <<"Physical Surface (\"Omega\") = {100};\n"
            ;

        std::ostringstream nameStr;
        nameStr.precision( 3 );
        nameStr << "StructuredGrid";
        gmshp->setPrefix( nameStr.str() );
        gmshp->setDescription( ostr.str() );
        return gmshp;
    }

    /*
     * \param filename : name of the file to be generated
     * \param outputs : vector containing outputs
     * \param parameter : vector containing parameter values
     * \param estimated_error : vector containing estimated error on outputs
     *
     */
    void generateGeoFileForOutputPlot(  vectorN_type outputs, vectorN_type parameter, vectorN_type estimated_error )
    {
        bool use_estimated_error=true;
        if( estimated_error(0) < 0 )
            use_estimated_error=false;

        Eigen::MatrixXf::Index index;
        double min_output = outputs.minCoeff(&index);
        //double min_scale=std::floor(min_output);
        double min_scale=min_output;
        double x=0;
        double output=0;
        double estimated_down=0;
        double estimated_up=0;
        double delta=0;

        std::string plotFile = "GMSH-outputs.geo";
        std::ofstream file_outputs_geo_gmsh ( plotFile, std::ios::out );
        file_outputs_geo_gmsh << "View \"outputs\" {\n";
        for(int i=0; i<outputs.size(); i++)
        {
            if( use_estimated_error )
            {
                delta=estimated_error(i);
                x=parameter(i);
                estimated_down=outputs(i);
                output=outputs(i)+delta/2.0;
                estimated_up=estimated_down+delta;
                file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
                file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<estimated_down<<", "<<min_scale<<"};\n";
                file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<estimated_up<<", "<<min_scale<<"};\n";
                file_outputs_geo_gmsh <<std::setprecision(14)<< "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
            }
            else
            {
                x=parameter(i);
                output=outputs(i);
                file_outputs_geo_gmsh << "SP("<<x<<",0,0){"<<output<<", "<<min_scale<<"};\n";
            }
        }

        std::string conclude=" }; \n ";
        conclude += "vid = PostProcessing.NbViews-1;\n";
        conclude += "View[vid].Axes = 1;\n";
        conclude += "View[vid].Type = 2;\n\n";
        conclude += "For i In {0:vid-1}\n";
        conclude += "  View[i].Visible=0;\n";
        conclude += "EndFor\n";

        file_outputs_geo_gmsh<<conclude;
        file_outputs_geo_gmsh.close();

        /* Adds the generated file for automatic loading in Gmsh */
        Environment::olLoadInGmsh(plotFile);


    }


    /**
     * specific interface for OpenTURNS
     *
     * \param X input vector of size N
     * \param N size of input vector X
     * \param Y input vector of size P
     * \param P size of input vector Y
     */
    virtual void run( const double * X, unsigned long N, double * Y, unsigned long P )
    {
        if( this->worldComm().isMasterRank() )
        {
            std::cout<<"**************************************************"<<std::endl;
            std::cout<<"** You are using the function run whereas       **"<<std::endl;
            std::cout<<"** your model has not implemented this function **"<<std::endl;
            std::cout<<"**************************************************"<<std::endl;
        }
    }

    virtual value_type output( int output_index, parameter_type const& mu , element_type &u , bool need_to_solve=false) = 0;

public:

    /**
     * Set the finite element space to \p Vh and then build the reduced basis
     * space from the finite element space.
     */
    virtual void setFunctionSpaces( functionspace_ptrtype const& Vh )
    {
        Xh = Vh;
        XN->setFunctionSpace( Xh );
        //XN->setModel( this->shared_from_this() );
    }

    /**
     * \brief Returns the function space
     */
    mesh_ptrtype const& mesh() const
    {
        return Xh->mesh();
    }

    /**
     * \brief Returns the function space
     */
    space_ptrtype const& functionSpace() const
    {
        return Xh;
    }

    virtual
    typename functionspace_type::mesh_support_vector_type
    functionspaceMeshSupport( mesh_ptrtype const& mesh ) const
        {
            return typename functionspace_type::mesh_support_vector_type();
        }

    /**
     * \brief Returns the reduced basis function space
     */
    rbfunctionspace_ptrtype const& rBFunctionSpace() const
    {
        return XN;
    }


    //! return the parameter space
    parameterspace_ptrtype const& parameterSpace() const
    {
        return Dmu;
    }


    //! Return a P1 vector field to warp the domain if we want to visualize the
    //! the effect of geometric parameters.
    virtual displacement_field_ptrtype meshDisplacementField( parameter_type const& mu )
    {
        CHECK( false ) << "Error : the function fieldForWarp() has been called without beeing rewrite.\n";

        displacement_field_ptrtype u;
        return u;
    }

    bool hasDisplacementField() const { return M_has_displacement_field; }
    void setHasDisplacementField( bool const _hwf ) { M_has_displacement_field=_hwf; }


protected :
    std::string M_prefix;

    parameterspace_ptrtype Dmu;
    functionspace_ptrtype Xh;
    rbfunctionspace_ptrtype XN;

    backend_ptrtype M_backend;

    CRBModelDB M_crbModelDb;

    funs_type M_funs;
    funsd_type M_funs_d;

    deim_vector_type M_deims;
    mdeim_vector_type M_mdeims;
    bool M_is_initialized, M_has_displacement_field;

    sparse_matrix_ptrtype M;
    sparse_matrix_ptrtype M_energy_matrix;
    sparse_matrix_ptrtype M_mass_matrix;

    operatorcomposite_ptrtype M_compositeA;
    operatorcomposite_ptrtype M_compositeM;
    std::vector< functionalcomposite_ptrtype > M_compositeF;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<sparse_matrix_ptrtype> M_Mq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Jqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_linearAqm;
    std::vector< std::vector<sparse_matrix_ptrtype> > M_Mqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;
    std::vector< std::vector<std::vector<vector_ptrtype> > > M_Rqm;

    sparse_matrix_ptrtype M_monoA;
    sparse_matrix_ptrtype M_monoM;
    std::vector<vector_ptrtype> M_monoF;

    beta_vector_light_type M_betaAq;
    beta_vector_light_type M_betaMq;
    std::vector<beta_vector_light_type> M_betaFq;

    std::vector< std::string > M_betaAqString;
    std::vector< std::string > M_betaMqString;
    std::vector< std::vector< std::string > > M_betaFqString;
    std::vector< Expr<GinacEx<2> > > M_ginacAq;
    std::vector< Expr<GinacEx<2> > > M_ginacMq;
    std::vector< std::vector< Expr<GinacEx<2> > > > M_ginacFq;
    std::vector< std::string > M_symbols_vec;
    std::string M_symbolicExpressionBuildDir;

    std::vector<double> M_betaTri;
    beta_vector_type M_betaAqm;
    beta_vector_type M_betaMqm;
    std::vector<beta_vector_type> M_betaFqm;

    beta_vector_type M_betaJqm;
    std::vector<beta_vector_type> M_betaRqm;

    beta_vector_type M_betaInitialGuess;

    lhs_light_type M_lhs;
    lhs_light_type M_mass;
    rhs_light_type M_rhs;
    rhs_light_type M_output;

    std::map<int,double> M_eim_error_mq;
    std::map<int,double> M_eim_error_aq;
    std::vector< std::map<int,double> > M_eim_error_fq;

    bool M_rebuildDb;
    int M_dbLoad;
    std::string M_dbFilename;
    std::string M_dbId;
    int M_dbUpdate;
    int M_serEimFreq;
    int M_serRbFreq;
    bool M_serErrorEstimation;
    bool M_crbUseNewton;

    bool M_isOnlineModel;
private :
    std::map<std::string,std::string > M_additionalModelFiles;
};




}//Feel
#endif /* __Model_H */
