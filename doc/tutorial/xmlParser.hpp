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

using namespace std;

#define MY_ENCODING "UTF-8"

#define OUTPUT_ATTRIBUTE 0
#define DISCRETE_ATTRIBUTE 1
#define CONTINUOUS_ATTRIBUTE 2

#define STR(a) (new string(a))

namespace xmlParse {

	class parameter {
		private:
			string* _name;
			int _type;
			string* _cmdName;
			string* _latex;
			string* _values;

		public:
			parameter(string* name, int type, string* cmdName, string* latex, string* values) {
				_name = name;
				_type = type;
				if (_cmdName!=NULL)
					_cmdName=cmdName;
				else
					_cmdName=name;
				_latex=latex;
				_values = values;
			}
			inline vector<string> getAttrNames() {
				vector<string> attrNames;
				attrNames.push_back("name");
				if (_type!=0)
					attrNames.push_back("type");
				if (_cmdName!=0)
					attrNames.push_back("cmd_name");
				if (_latex!=0)
					attrNames.push_back("latex");
				return attrNames;
			}
			inline vector<string> getAttrValues() {
				vector<string> attrValues;
				attrValues.push_back(*_name);
				if (_type!=0)
					attrValues.push_back( ( _type == 1 ? "discrete" : "continuous" ) );
				if (_cmdName!=0)
					attrValues.push_back(*_cmdName);
				if (_latex!=0)
					attrValues.push_back(*_latex);
				return attrValues;
			}
			inline string getValues() {
				if (_values!=0)
					return *_values;
				else
					return "";
			}
			inline string getName() {
				return *_name;
			}
	};

	class output : public parameter {
		private:
			vector<parameter> _dependencies;
			vector<string> _funcs;

		public:
			output(string* name, string* latex, vector<parameter> dependencies, vector<string> funcs) : parameter(name, 0, NULL, latex, NULL) {
				_dependencies=dependencies;
				_funcs=funcs;
			}
			inline vector<parameter> getDependencies() {
				return _dependencies;
			}
			inline vector<string> getFuncs() {
				return _funcs;
			}
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
					xmlNodePtr paramNode = createNode("output",outputs[i].getAttrNames(),outputs[i].getAttrValues(),outputs[i].getValues());
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
