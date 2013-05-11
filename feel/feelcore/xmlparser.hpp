/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


  This file is part of the Feel library

  Author(s):
  Florent Vielfaure <florent.vielfaure@gmail.com>

  Date: 2009-02-22

  Copyright (C) 2009 Université de Grenoble 1

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
   \file xmlParser.hpp
   \author Florent Vielfaure <florent.vielfaure@gmail.com>
   \date 2009-02-22
 */
#if !defined(XML_PARSER_H)
#define XML_PARSER_H 1

#include <vector>
#include <cstdio>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include <boost/parameter.hpp>
#include <feel/feelcore/parameter.hpp>

#define MY_ENCODING "UTF-8"

#define OUT_ATTR 0
#define DISC_ATTR 1
#define CONT_ATTR 2

namespace Feel
{

class Parameter_impl
{
private:
    std::string* name;
    int type;
    std::string* cmdName;
    std::string* latex;
    std::string* values;

public:
    Parameter_impl() {}
    template <class ArgumentPack>
    Parameter_impl( ArgumentPack const& args )
    {
        name = new std::string( ( const char* )args[_name] );
        type = ( int )args[_type];
        cmdName = ( args[_cmdName | 0] ?
                    new std::string( ( const char* )args[_cmdName | ( const char* )0] ) :
                    new std::string( *name ) );
        latex = ( args[_latex | 0] ? new std::string( ( const char* )args[_latex | ( const char* )0] ) : ( std::string* )0 );
        values = ( args[_values | 0] ? new std::string( ( const char* )args[_values | ( const char* )0] ) : ( std::string* )0 );
    }
    /**
     * \brief Function used to know the attributes available in the parameter
     * \return the name of the attributes that are given for the parameter
     */
    inline std::vector<std::string> getAttrNames()
    {
        std::vector<std::string> attrNames;
        attrNames.push_back( "name" );

        if ( type )
            attrNames.push_back( "type" );

        if ( cmdName )
            attrNames.push_back( "cmd_name" );

        if ( latex )
            attrNames.push_back( "latex" );

        return attrNames;
    }
    /**
     * \brief Used to retrieve the attributes values for the parameter.
     * The array contain the values of the parameters in the same order than the
     * attribute names returned in getAttrNames() function.
     * \return the values of the attributes for a parameter
     */
    inline std::vector<std::string> getAttrValues()
    {
        std::vector<std::string> attrValues;
        attrValues.push_back( *name );

        if ( type )
            attrValues.push_back( ( type == 1 ? "discrete" : "continuous" ) );

        if ( cmdName )
            attrValues.push_back( *cmdName );

        if ( latex )
            attrValues.push_back( *latex );

        return attrValues;
    }
    /**
     * \return the values of the parameter
     */
    inline std::string getValues()
    {
        return *values;
    }
    /**
     * \return the name of the parameter
     */
    inline std::string getName()
    {
        return *name;
    }
};

/**
 * \class Parameter
 * \brief parameter class to describe code inputs
 */
class Parameter : public Parameter_impl
{
public:
    Parameter() {}
    BOOST_PARAMETER_CONSTRUCTOR(
        Parameter, ( Parameter_impl ), tag
        , ( required ( name,* ) ) ( required ( type,* ) ) ( optional ( cmdName,* ) ) ( optional ( latex,* ) ) ( optional ( values,* ) ) )
};

class Output_impl : public Parameter
{
private:
    std::vector<Parameter> dependencies;
    std::vector<std::string> funcs;

public:
    template <class ArgumentPack>
    Output_impl( ArgumentPack const& args )
        :
        Parameter( _name=args[_name],
                   _type=0,
                   _latex=args[_latex ] )
    {
        dependencies=args[_dependencies | std::vector<Parameter>() ];
        funcs=args[_funcs | std::vector<std::string>() ];
    }
    /**
     * \return the dependencies of the output
     */
    inline std::vector<Parameter> getDependencies()
    {
        return dependencies;
    }
    /**
     * \return the functions associated to a dependency of an output
     */
    inline std::vector<std::string> getFuncs()
    {
        return funcs;
    }
};

/**
 * \class Output
 * \brief output class to describe code outputs
 */
class Output : public Output_impl
{
public:
    BOOST_PARAMETER_CONSTRUCTOR(
        Output,
        ( Output_impl ),
        tag,
        ( required ( name,* ) )
        ( optional ( latex,* ) )
        ( optional ( dependencies,* ) )
        ( optional ( funcs,* ) )
    )
};

class xmlParser
{
public:
    /**
     * writeResponse
     *
     * Writes the options available in the c++ code.
     **/
    static void writeResponse( std::string filename,std::string name,std::vector<Parameter> params,std::vector<Output> Outputs )
    {
        xmlDocPtr doc = xmlNewDoc( ( xmlChar* ) "1.0" );
        xmlNodePtr rootNode = xmlNewNode( 0, ( xmlChar* ) "response" );
        xmlDocSetRootElement( doc, rootNode );
        xmlNodePtr progNode = createNode( "program_name",name );
        xmlAddChild( rootNode,progNode );

        for ( unsigned int i=0; i<params.size(); ++i )
        {
            xmlNodePtr paramNode = createNode( "param",params[i].getAttrNames(),params[i].getAttrValues(),params[i].getValues() );
            xmlAddChild( rootNode,paramNode );
        }

        for ( unsigned int i=0; i<Outputs.size(); ++i )
        {
            xmlNodePtr paramNode = createNode( "output",Outputs[i].getAttrNames(),Outputs[i].getAttrValues() );

            for ( unsigned int j=0; j<Outputs[i].getDependencies().size(); ++j )
            {
                xmlNode* aNewNode = xmlNewNode( 0,( xmlChar* ) "depend" );
                xmlSetProp( aNewNode, ( xmlChar* ) "value", ( xmlChar* ) Outputs[i].getDependencies()[j].getName().c_str() );
                xmlNodeSetContent( aNewNode,( xmlChar* ) Outputs[i].getFuncs()[j].c_str() );
                xmlAddChild( paramNode,aNewNode );
            }

            xmlAddChild( rootNode,paramNode );
        }

        xmlSaveFileEnc( filename.c_str(), doc, MY_ENCODING );

        xmlFreeDoc( doc );
    }

    /**
     * writeResult
     *
     * Appends the result in the xml file.
     *
     **/
    static void writeResult( std::string filename,
                             std::string name,
                             std::vector<Parameter> params,
                             std::vector<Output> Outputs,
                             std::vector<std::string> paramValues,
                             std::vector<std::string> OutputValues )
    {
        xmlDoc *doc = 0;
        xmlNode *root_element = 0;

        doc = xmlReadFile( filename.c_str(), 0, 0 );

        if ( doc == 0 )
        {
            doc = xmlNewDoc( ( xmlChar* ) "1.0" );
            root_element = xmlNewNode( 0, ( xmlChar* ) "result" );
            xmlDocSetRootElement( doc, root_element );
            printf( "[xmlParser] warning : the file %s has been generated\n", filename.c_str() );
        }

        root_element = xmlDocGetRootElement( doc );

        if ( root_element == 0 )
        {
            printf( "[xmlParser] error : the file %s is empty\n", filename.c_str() );
            return;
        }

        xmlNode* aNode=findNode( root_element,"program", name.c_str() );
        xmlNode* aNewNode;

        for ( unsigned int i=0; i<params.size(); ++i )
        {
            aNewNode=findNode( aNode, params[i].getName(), paramValues[i] );
            aNode=aNewNode;
        }

        for ( unsigned int i=0; i<Outputs.size(); ++i )
        {
            xmlNode *cur_node=0;

            for ( cur_node = aNode->children; cur_node; cur_node = cur_node->next )
            {
                if ( ( cur_node->type == XML_ELEMENT_NODE ) &&
                        ( std::strcmp( ( char* )cur_node->name,Outputs[i].getName().c_str() )==0 ) )
                {
                    xmlUnlinkNode( cur_node );
                    xmlFreeNode( cur_node );
                }
            }

            aNewNode=createNode( Outputs[i].getName(), OutputValues[i] );
            xmlAddChild( aNode,aNewNode );
        }

        xmlSaveFileEnc( filename.c_str(), doc, MY_ENCODING );
    }
private:
    static xmlNode* findNode( xmlNode* aNode, std::string nodeName, std::string attrVal )
    {
        xmlNode *cur_node = 0;

        for ( cur_node = aNode->children; cur_node; cur_node = cur_node->next )
        {
            if ( ( cur_node->type == XML_ELEMENT_NODE ) &&
                    ( std::strcmp( ( char* )cur_node->name,nodeName.c_str() )==0 ) )
            {
                xmlChar* val=xmlGetProp( cur_node,( xmlChar* ) "value" );

                if ( ( val!=0 ) && ( std::strcmp( ( char* )val,attrVal.c_str() )==0 ) )
                {
                    return cur_node;
                }
            }
        }

        xmlNode* aNewNode = xmlNewNode( 0,( xmlChar* ) nodeName.c_str() );
        xmlSetProp( aNewNode, ( xmlChar* ) "value", ( xmlChar* ) attrVal.c_str() );
        xmlAddChild( aNode,aNewNode );
        return aNewNode;
    }

    static xmlNode* createNode( std::string nodeName, std::vector<std::string> attrNames, std::vector<std::string> attrValues, std::string value="" )
    {
        if ( attrNames.size()!=attrValues.size() )
        {
            printf( "[xmlParser] error : attributes names list size != attributes values list size\n" );
            return 0;
        }

        xmlNode* aNewNode = xmlNewNode( 0,( xmlChar* ) nodeName.c_str() );

        for ( unsigned int i=0; i<attrNames.size(); ++i )
            xmlSetProp( aNewNode, ( xmlChar* ) attrNames[i].c_str(), ( xmlChar* ) attrValues[i].c_str() );

        if ( std::strcmp( value.c_str(),"" )!=0 )
            xmlNodeSetContent( aNewNode,( xmlChar* ) value.c_str() );

        return aNewNode;
    }

    static xmlNode* createNode( std::string nodeName, std::string value="" )
    {
        xmlNode* aNewNode = xmlNewNode( 0,( xmlChar* ) nodeName.c_str() );

        if ( std::strcmp( value.c_str(),"" )!=0 )
            xmlNodeSetContent( aNewNode,( xmlChar* ) value.c_str() );

        return aNewNode;
    }
};
}

#endif /* XML_PARSER_H */
