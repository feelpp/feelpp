/* -*- mode: c++ -*-


  This file is part of the Life library

  Author(s):
  Florent Vielfaure <florent.vielfaure@gmail.com>

  Date: 2009-02-22

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

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

using namespace std;

#define MY_ENCODING "UTF-8"

#define OUT_ATTR 0
#define DISC_ATTR 1
#define CONT_ATTR 2

namespace xmlParse {

BOOST_PARAMETER_NAME(name)
BOOST_PARAMETER_NAME(type)
BOOST_PARAMETER_NAME(latex)
BOOST_PARAMETER_NAME(cmdName)
BOOST_PARAMETER_NAME(values)
BOOST_PARAMETER_NAME(dependencies)
BOOST_PARAMETER_NAME(funcs)

class parameter_impl {
private:
    string* name;
    int type;
    string* cmdName;
    string* latex;
    string* values;

public:
    parameter_impl() {}
    template <class ArgumentPack>
    parameter_impl(ArgumentPack const& args) {
        name = new string((const char*)args[_name]);
        type = (int)args[_type];
        cmdName = ( args[_cmdName | NULL] ?
                    new string((const char*)args[_cmdName | NULL]) :
                    new string(*name) );
        latex = (args[_latex | NULL] ? new string((const char*)args[_latex | NULL]) : NULL);
        values = (args[_values | NULL] ? new string((const char*)args[_values | NULL]) : NULL);

        /* printf("name=%s\n",name->c_str());
        printf("type=%d\n",type);
        printf("cmdName=%s\n",cmdName->c_str());
        if (latex)
            printf("latex=%s\n",latex->c_str());
        if (values)
        printf("values=%s\n",values->c_str()); */
    }
    inline vector<string> getAttrNames() {
        // printf("getAttrNames()\n");
        vector<string> attrNames;
        attrNames.push_back("name");
        if (type)
            attrNames.push_back("type");
        if (cmdName)
            attrNames.push_back("cmd_name");
        if (latex)
            attrNames.push_back("latex");
        return attrNames;
    }
    inline vector<string> getAttrValues() {
        // printf("getAttrValues()\n");
        vector<string> attrValues;
        attrValues.push_back(*name);
        if (type)
            attrValues.push_back( ( type == 1 ? "discrete" : "continuous" ) );
        if (cmdName)
            attrValues.push_back(*cmdName);
        if (latex)
            attrValues.push_back(*latex);
        return attrValues;
    }
    inline string getValues() {
        // printf("getValues()\n");
        return *values;
    }
    inline string getName() {
        // printf("getName()\n");
        return *name;
    }
};

class parameter : public parameter_impl {
    public:
    parameter() {}
    BOOST_PARAMETER_CONSTRUCTOR(
        parameter, (parameter_impl), tag
        , (required (name,*)) (required (type,*)) (optional (cmdName,*)) (optional (latex,*)) (optional (values,*)))
};

class output_impl : public parameter {
private:
    vector<parameter> dependencies;
    vector<string> funcs;

public:
    template <class ArgumentPack>
    output_impl(ArgumentPack const& args) : parameter(_name=args[_name], _type=0, _latex=args[_latex]) {
        dependencies=args[_dependencies];
        funcs=args[_funcs];
    }
    inline vector<parameter> getDependencies() {
        return dependencies;
    }
    inline vector<string> getFuncs() {
        return funcs;
    }
};

class output : public output_impl {
    public:
    BOOST_PARAMETER_CONSTRUCTOR(
        output, (output_impl), tag
        , (required (name,*)) (optional (latex,*)) (required (dependencies,*)) (required (funcs,*)))
};

class xmlParser {
public:
    /**
     * writeResponse
     *
     * Writes the options available in the c++ code.
     **/
    static void writeResponse(string filename,string name,vector<parameter> params,vector<output> outputs) {
        xmlDocPtr doc = xmlNewDoc((xmlChar*) "1.0");
        xmlNodePtr rootNode = xmlNewNode(NULL, (xmlChar*) "response");
        xmlDocSetRootElement(doc, rootNode);
        xmlNodePtr progNode = createNode("program_name",name);
        xmlAddChild(rootNode,progNode);

        for (unsigned int i=0; i<params.size(); ++i) {
            xmlNodePtr paramNode = createNode("param",params[i].getAttrNames(),params[i].getAttrValues(),params[i].getValues());
            xmlAddChild(rootNode,paramNode);
        }
        for (unsigned int i=0; i<outputs.size(); ++i) {
            xmlNodePtr paramNode = createNode("output",outputs[i].getAttrNames(),outputs[i].getAttrValues());
            for (unsigned int j=0; j<outputs[i].getDependencies().size(); ++j) {
                xmlNode* aNewNode = xmlNewNode(NULL,(xmlChar*) "depend");
                xmlSetProp(aNewNode, (xmlChar*) "value", (xmlChar*) outputs[i].getDependencies()[j].getName().c_str());
                xmlNodeSetContent(aNewNode,(xmlChar*) outputs[i].getFuncs()[j].c_str());
                xmlAddChild(paramNode,aNewNode);
            }
            xmlAddChild(rootNode,paramNode);
        }
        xmlSaveFileEnc(filename.c_str(), doc, MY_ENCODING);

        xmlFreeDoc(doc);
    }

    /**
     * writeResult
     *
     * Appends the result in the xml file.
     *
     **/
    static void writeResult(string filename,
                            string name,
                            vector<parameter> params,
                            vector<output> outputs,
                            vector<string> paramValues,
                            vector<string> outputValues) {
        xmlDoc *doc = NULL;
        xmlNode *root_element = NULL;

        doc = xmlReadFile(filename.c_str(), NULL, 0);

        if (doc == NULL) {
            doc = xmlNewDoc((xmlChar*) "1.0");
            root_element = xmlNewNode(NULL, (xmlChar*) "result");
            xmlDocSetRootElement(doc, root_element);
            printf("[xmlParser] warning : the file %s has been generated\n", filename.c_str());
        }

        root_element = xmlDocGetRootElement(doc);

        if (root_element == NULL) {
            printf("[xmlParser] error : the file %s is empty\n", filename.c_str());
            return;
        }

        xmlNode* aNode=findNode(root_element,"program", name.c_str());
        xmlNode* aNewNode;
        for (unsigned int i=0; i<params.size();++i) {
            aNewNode=findNode(aNode, params[i].getName(), paramValues[i]);
            aNode=aNewNode;
        }
        for (unsigned int i=0; i<outputs.size();++i) {
            xmlNode *cur_node=NULL;
            for (cur_node = aNode->children; cur_node; cur_node = cur_node->next) {
                if ( (cur_node->type == XML_ELEMENT_NODE) &&
                     (std::strcmp((char*)cur_node->name,outputs[i].getName().c_str())==0) ) {
                    xmlUnlinkNode(cur_node);
                    xmlFreeNode(cur_node);
                }
            }
            aNewNode=createNode(outputs[i].getName(), outputValues[i]);
            xmlAddChild(aNode,aNewNode);
        }

        xmlSaveFileEnc(filename.c_str(), doc, MY_ENCODING);
    }
private:
    static xmlNode* findNode(xmlNode* aNode, string nodeName, string attrVal) {
        xmlNode *cur_node = NULL;

        for (cur_node = aNode->children; cur_node; cur_node = cur_node->next) {
            if ( (cur_node->type == XML_ELEMENT_NODE) &&
                 (std::strcmp((char*)cur_node->name,nodeName.c_str())==0) ) {
                xmlChar* val=xmlGetProp(cur_node,(xmlChar*) "value");
                if ( (val!=NULL) && (std::strcmp((char*)val,attrVal.c_str())==0) ) {
                    return cur_node;
                }
            }
        }
        xmlNode* aNewNode = xmlNewNode(NULL,(xmlChar*) nodeName.c_str());
        xmlSetProp(aNewNode, (xmlChar*) "value", (xmlChar*) attrVal.c_str());
        xmlAddChild(aNode,aNewNode);
        return aNewNode;
    }

    static xmlNode* createNode(string nodeName, vector<string> attrNames, vector<string> attrValues, string value="") {
        if (attrNames.size()!=attrValues.size()) {
            printf("[xmlParser] error : attributes names list size != attributes values list size\n");
            return NULL;
        }
        xmlNode* aNewNode = xmlNewNode(NULL,(xmlChar*) nodeName.c_str());
        for (unsigned int i=0; i<attrNames.size(); ++i)
            xmlSetProp(aNewNode, (xmlChar*) attrNames[i].c_str(), (xmlChar*) attrValues[i].c_str());
        if (std::strcmp(value.c_str(),"")!=0)
            xmlNodeSetContent(aNewNode,(xmlChar*) value.c_str());
        return aNewNode;
    }

    static xmlNode* createNode(string nodeName, string value="") {
        xmlNode* aNewNode = xmlNewNode(NULL,(xmlChar*) nodeName.c_str());
        if (std::strcmp(value.c_str(),"")!=0)
            xmlNodeSetContent(aNewNode,(xmlChar*) value.c_str());
        return aNewNode;
    }
	};
}

#endif /* XML_PARSER_H */
