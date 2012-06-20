/*
 *  codeare Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301  USA
 */

#ifndef __CONFIGURABLE_HPP__
#define __CONFIGURABLE_HPP__

#include "tinyxml.h"
#include "tinyxml/xpath_static.h"

#include <string>

using namespace TinyXPath;

/**
 * @brief Skeleton of an XML configurable class
 */
class Configurable {



 public:


	/**
	 * @brief           Construct and initialise configuration DOM
	 */
	Configurable() {

		m_config_doc             = new TiXmlDocument    ();
		m_config_doc->LinkEndChild(new TiXmlDeclaration ("1.0", "", ""));
		m_config_doc->LinkEndChild(new TiXmlElement     ("config"));
		
	}


	/**
	 * @brief           Delete configuration DOM
	 */
	virtual ~Configurable() {
		delete m_config_doc;
	}

	
	/**
	 * @brief           Set a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, const char* value) {
		Configuration()->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, const std::string value) {
		Configuration()->SetAttribute (name, value.c_str());
	}

	
	/**
	 * @brief           Set a integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, int value) {
		Configuration()->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a bool type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, bool value) {
		Configuration()->SetAttribute (name, (int)value);
	}

	
	/**
	 * @brief           Set a size type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, size_t value) {
		Configuration()->SetAttribute (name, (int)value);
	}

	
	/**
	 * @brief           Set a float type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void 
	SetAttribute        (const char* name, double value) {
		Configuration()->SetDoubleAttribute (name, value);
	}


	/**
	 * @brief           Get a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @return          String representation of value
	 */
	inline const char*
	Attribute           (const char* name) const {
		return Configuration()->Attribute (name);
	}
	

	/**
	 * @brief           Get a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Length
	 */
	inline int
	Attribute           (const char* name, std::string* value) const {
		value = new std::string (Configuration()->Attribute (name));
		return value->length();
	}
	

	/**
	 * @brief           Get an integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, int* value) const {
		return Configuration()->QueryIntAttribute (name, value);
	}

	
	/**
	 * @brief           Set a bool type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, size_t* value) const {
		int ival, success;
		success = Configuration()->QueryIntAttribute (name, &ival);
		*value = (size_t) ival;
		return success;
	}

	
	/**
	 * @brief           Set a bool type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, bool* value) const {
		int ival, success;
		success = Configuration()->QueryIntAttribute (name, &ival);
		*value = (bool) ival;
		return success;
	}

	
	/**
	 * @brief           Set a double type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, double* value) const {
		return Configuration()->QueryDoubleAttribute (name, value);
	}


	/**
	 * @brief           Set a float type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, float* value) const {
		return Configuration()->QueryDoubleAttribute (name, (double*)value);
	}


	/**
	 * @brief           Get text of a node
	 * 
	 * @param  path     Path to node
	 * @return          Value
	 */
	inline const char* 
	GetText             (const char* path) const {
		return ((TiXmlElement*) TinyXPath::XNp_xpath_node (Configuration(), path))->GetText();
	}

		
	/**
	 * @brief           Get a particular node
	 * 
	 * @param  path     Path to node
	 * @return          Value
	 */
	inline TiXmlElement*
	GetElement          (const char* path) const {
		return (TiXmlElement*) TinyXPath::XNp_xpath_node (Configuration(), path);
	}

		
	/**
	 * @brief Serialize configuration to string
	 *
	 * @return          String representation of configuration
	 */
	inline const char*
	GetConfig           ()       const {

		TiXmlPrinter printer;
		std::string  tstr;
		char*        tcstr;

		m_config_doc->Accept( &printer );
		tstr  = printer.CStr();
		tcstr = new char[tstr.length() + 1];
		strcpy (tcstr, tstr.c_str());

		return tcstr;

	}


	/**
	 * @brief          Set configuration from string
	 *
	 * @param  cstr    Serialized configuration
	 */
	inline void 
	SetConfig          (const char* cstr) {

		m_config_doc = new TiXmlDocument();
		m_config_doc->Parse(cstr);

	}


	/**
	 * @brief          Get configuration element for manipulation. (lhs)
	 *
	 * @return         Pointer root element "config"
	 */
	TiXmlElement* 
	Configuration      () {
		return m_config_doc->RootElement();
	}


	/**
	 * @brief          Get configuration element for reading (rhs)
	 *
	 * @return         Const pointer to Element "Config"
	 */
	const TiXmlElement* 
	Configuration      () const {
		return m_config_doc->RootElement();
	}


	/**
	 * @brief          Dump XML configuration to file
	 *
	 * @param  fname   Name of output file
	 */
	inline void 
	DumpConfig        (const char* fname) {
		m_config_doc->SaveFile (fname);
	}
	

	/**
	 * @brief          Read configuration from XML file
	 *
	 * @param  fname   Name of input file
	 */
	void 
	ReadConfig        (const char* fname) {
		m_config_doc->LoadFile (fname);
	}
	

	/**
	 * @brief          Read configuration from XML file
	 *
	 * @param  file    File pointer
	 */
	void 
	ReadConfig        (FILE* file) {
		m_config_doc->LoadFile (file);
	}
	


 protected:
	
	TiXmlDocument*      m_config_doc;  /**< @brief DOM docuemnt holding the configuration */ 

};



#endif
