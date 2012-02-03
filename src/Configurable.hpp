/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
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
#include <string>

/**
 * @brief Skeleton of an XML configurable class
 */
class Configurable {



 public:


	/**
	 * @brief           Construct and initialise configuration DOM
	 */
	Configurable() {

		std::cout << "We are never called?" << std::endl;

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
		m_config_doc->RootElement()->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, const std::string value) {
		m_config_doc->RootElement()->SetAttribute (name, value.c_str());
	}

	
	/**
	 * @brief           Set a integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, int value) {
		m_config_doc->RootElement()->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a bool type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, bool value) {
		m_config_doc->RootElement()->SetAttribute (name, (int)value);
	}

	
	/**
	 * @brief           Set a size type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute        (const char* name, size_t value) {
		m_config_doc->RootElement()->SetAttribute (name, (int)value);
	}

	
	/**
	 * @brief           Set a float type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void 
	SetAttribute        (const char* name, double value) {
		m_config_doc->RootElement()->SetDoubleAttribute (name, value);
	}


	/**
	 * @brief           Get a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @return          String representation of value
	 */
	inline const char*
	Attribute           (const char* name) const {
		return m_config_doc->RootElement()->Attribute (name);
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
		value = new std::string (m_config_doc->RootElement()->Attribute (name));
		return value->length();
	}
	

	/**
	 * @brief           Set an integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 * @return          Status
	 */
	inline int
	Attribute           (const char* name, int* value) const {
		return m_config_doc->RootElement()->QueryIntAttribute (name, value);
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
		return m_config_doc->RootElement()->QueryIntAttribute (name, (int*)value);
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
		return m_config_doc->RootElement()->QueryIntAttribute (name, (int*)value);
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
		return m_config_doc->RootElement()->QueryDoubleAttribute (name, value);
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
		return m_config_doc->RootElement()->QueryDoubleAttribute (name, (double*)value);
	}


	/**
	 * @brief           Serialize configuration to string
	 *
	 * @return          String representation of configuration
	 */
	inline const char*
	GetConfig           ()             {

		int err = 0;
		TiXmlPrinter printer;
		m_config_doc->Accept( &printer );
		std::string temp  = printer.CStr();
		char* t = new char[temp.length() + 1];
		strcpy (t, temp.c_str());
		return t;

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
	 * @brief          Get configuration element for more elaborate manipulation. (i.e. Adding children etc)
	 *
	 * @return         Reference to Element "Config"
	 */
	TiXmlElement* 
	Configuration      () {
		return m_config_doc->RootElement();
	}


	/**
	 * @brief          Dump XML configuration to file
	 *
	 * @param  fname   Name of output file
	 */
	void 
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
	
	TiXmlDocument*      m_config_doc;  /**< Containing XML document              */ 

};



#endif
