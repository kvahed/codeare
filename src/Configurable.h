#ifndef __CONFIGURABLE_H__
#define __CONFIGURABLE_H__

#include "tinyxml.h"


/**
 * @brief Skeleton of an XML configurable class
 */
class Configurable {



 public:


	/**
	 * @brief           Initialise configuration 
	 */
	Configurable() {

		m_config_doc  = new TiXmlDocument();
		m_config_decl = new TiXmlDeclaration ("1.0", "", "");
		m_config      = new TiXmlElement     ("Config");

		m_config_doc->LinkEndChild( m_config_decl );
		m_config_doc->LinkEndChild( m_config );
		
	}


	/**
	 * @brief           Delete configuration
	 */
	virtual ~Configurable() {
		delete m_config_doc;
	};

	
	/**
	 * @brief           Set a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute           (const char* name, const char* value) {
		m_config->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void
	SetAttribute           (const char* name, int value) {
		m_config->SetAttribute (name, value);
	}

	
	/**
	 * @brief           Set a float type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline void 
	SetAttribute           (const char* name, double value) {
		m_config->SetDoubleAttribute (name, value);
	}


	/**
	 * @brief           Set a string type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline const char*
	Attribute           (const char* name) const {
		return m_config->Attribute (name);
	}
	

	/**
	 * @brief           Set a integer type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline const char*
	Attribute           (const char* name, int* value) const {
		return m_config->Attribute (name, value);
	}

	
	/**
	 * @brief           Set a float type attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	inline const char*
	Attribute           (const char* name, double* value) const {
		return m_config->Attribute (name, value);
	}


	/**
	 * @brief           Serialize configuration to string
	 */
	inline void 
	GetConfig           (char* cstr) {

		std::string temp = "";
		temp << *(m_config_doc);
		strcpy (cstr, temp.c_str());

	}


	/**
	 * @brief           Serialize configuration to string
	 */
	inline const char*
	GetConfig           ()             {

		std::string temp = "";
		temp << *(m_config_doc);
		char* t = new char[temp.length()];
		strcpy (t, temp.c_str());
		return t;

	}


	/**
	 * @brief          Set configuration from string
	 */
	inline void 
	SetConfig          (const char* cstr) {

		m_config_doc->Clear();
		m_config_doc->Parse(cstr);
		m_config = m_config_doc->RootElement();

	}


	/**
	 * @brief          Get configuration element
	 */
	TiXmlElement* 
	Configuration      () {
		
		return m_config;
		
	}

	void 
	DumpConfig        (const char* fname) {
		m_config_doc->SaveFile (fname);
	}
	

 protected:
	

	TiXmlElement*       m_config;      /**< Flat Configuration XML element.      */ 
	TiXmlDocument*      m_config_doc;  /**< Containing XML document              */ 
	TiXmlDeclaration*   m_config_decl; /**< Declaration                          */

};



#endif
