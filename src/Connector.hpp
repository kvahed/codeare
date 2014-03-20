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

#ifndef __CONNECTOR_HPP__
#define __CONNECTOR_HPP__

#include "Matrix.hpp"
#include "LocalConnector.hpp"
#include "RemoteConnector.hpp"

namespace RRClient {

enum ConType {LOCAL, REMOTE};

/**
 * @brief Connector skeleton
 *        Abstraction layer for local or remote access to reconstruction schemes
 */
class Connector {

public:

	/**
	 * @brief       Construct with service name and debug level.
	 *
	 * @param name  Service name
	 * @param debug Trace level
	 */
    Connector  (int args = 0, char** argv = 0,
    		const char* name = 0, const char* debug = 0) {
    	m_ct = (name==0) ? LOCAL : REMOTE;
        m_conn = (m_ct==LOCAL) ?
        		(Connection*) new LocalConnector (args, argv, name, debug) :
        		(Connection*) new RemoteConnector (args, argv, name, debug);
	}
	
	
	/**
	 * @brief      Close connection
	 */
	virtual 
	~Connector () {
		delete m_conn;
	}


	/**
	 * @brief           Request data procession on remote service
	 *
	 *                  @see RRStrategy::ReconStrategy::Process()
	 *
	 * @param  name     Recon method
	 * @return          Error code
	 */ 
	virtual inline codeare::error_code              
	Process             (const char* name) {
		codeare::error_code ec = codeare::OK;
		switch (m_ct) {
		case  LOCAL: ec = (codeare::error_code) ( (LocalConnector*) m_conn)->Process(name); break;
		case REMOTE: ec = (codeare::error_code) ((RemoteConnector*) m_conn)->Process(name); break;
		default: break;
		}
		return ec;
	}
	
	
	/**
	 * @brief           Prepare backend
	 *
	 *                  @see RRStrategy::ReconStrategy::Prepare()
	 *
	 * @param  name     Recon method
	 * @return          Error code
	 */ 
	virtual inline codeare::error_code              
	Prepare             (const char* name) {
		codeare::error_code ec = codeare::OK;
		switch (m_ct) {
		case  LOCAL: ec = (codeare::error_code) ( (LocalConnector*) m_conn)->Prepare(name); break;
		case REMOTE: ec = (codeare::error_code) ((RemoteConnector*) m_conn)->Prepare(name); break;
		default: break;
		}
		return ec;
	}
	
	
	/**
	 * @brief           Initialise remote service
	 *
	 *                  @see RRStrategy::ReconStrategy::Init()
	 *
	 * @param  name     Recon method
	 * @return          Error code
	 */ 
	virtual inline codeare::error_code              
	Init                (const char* name, const char* config) {
		codeare::error_code ec = codeare::OK;
		switch (m_ct) {
		case  LOCAL: ec = (codeare::error_code) ( (LocalConnector*) m_conn)->Init(name, config); break;
		case REMOTE: ec = (codeare::error_code) ((RemoteConnector*) m_conn)->Init(name, config); break;
		default: break;
		}
		return ec;
	}
	
	
	/**
	 * @brief           Finalise remote service
	 *
	 *                  @see RRStrategy::ReconStrategy::Finalise()
	 *
	 * @param  name     Recon method
	 * @return          Error error
	 */ 
	virtual inline codeare::error_code              
	Finalise            (const char* name) {
		codeare::error_code ec = codeare::OK;
		switch (m_ct) {
		case  LOCAL: ec = (codeare::error_code) ( (LocalConnector*) m_conn)->Finalise(name); break;
		case REMOTE: ec = (codeare::error_code) ((RemoteConnector*) m_conn)->Finalise(name); break;
		default: break;
		}
		return ec;
	}
	
	
	/**
	 * @brief           Transmit measurement data to remote service
	 *
	 *                  @see Database::SetMatrix
	 *
	 * @param  name     Name
	 * @param  m        Matrix
	 */
	template <class S> inline void 
	SetMatrix           (const std::string& name, Matrix<S>& m) const {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->SetMatrix(name, m); break;
		case REMOTE: ((RemoteConnector*) m_conn)->SetMatrix(name, m); break;
		default: break;
		}
	}
	
	
	/**
	 * @brief           Retrieve manipulated data from remote service
	 *
	 *                  @see Database::GetMatrix
	 *
	 * @param  name     Name
	 * @param  m        Receive storage
	 */
	template <class S> inline void 
	GetMatrix           (const std::string& name, Matrix<S>& m) const {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->GetMatrix(name, m); break;
		case REMOTE: ((RemoteConnector*) m_conn)->GetMatrix(name, m); break;
		default: break;
		}
	}
		
		
	/**
	 * @brief          Read configuration 
	 *
	 *                 @see Configurable::ReadConfig(const char* fname)
	 *                 @see Configurable::ReadConfig(FILE* file)
	 *
	 * @param config   Name of input file or file access pointer
	 */
	template <class S> inline void 
	ReadConfig        (S config) {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->ReadConfig(config); break;
		case REMOTE: ((RemoteConnector*) m_conn)->ReadConfig(config); break;
		default: break;
		}
	}
	

	/**
	 * @brief          Read configuration
	 *
	 *                 @see Configurable::ReadConfig(const char* fname)
	 *                 @see Configurable::ReadConfig(FILE* file)
	 *
	 * @param config   Name of input file or file access pointer
	 */
	void SetConfig        (const char *confstr) {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->SetConfig(confstr); break;
		case REMOTE: ((RemoteConnector*) m_conn)->SetConfig(confstr); break;
		default: break;
		}
	}


	/**
	 * @brief           Set a string type attribute
	 *
	 *                  @see Configurable::SetAttribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	template <class S> inline void
	SetAttribute        (const char* name, S value) {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->SetAttribute (name, value); break;
		case REMOTE: ((RemoteConnector*) m_conn)->SetAttribute (name, value); break;
		default: break;
		}
	}


	/**
	 * @brief           Set a string type attribute
	 *
	 *                  @see Configurable::Attribute
	 *
	 * @param  name     Attribute name 
	 * @param  value    Attribute value
	 */
	template <class S> inline int
	Attribute        (const char* name, S* value) {
		switch (m_ct) {
		case  LOCAL: ( (LocalConnector*) m_conn)->Attribute (name, value); break;
		case REMOTE: ((RemoteConnector*) m_conn)->Attribute (name, value); break;
		default: break;
		}
	}

	
	/**
	 * @brief           Get a string type attribute
	 *
	 *                  @see Configurable::Attribute
	 *
	 * @param  name     Attribute name 
	 * @return          Attribute value
	 */
	inline const char*
	Attribute          (const char* name) {
		const char* ret;
		switch (m_ct) {
		case  LOCAL: ret = ( (LocalConnector*) m_conn)->Attribute (name); break;
		case REMOTE: ret = ((RemoteConnector*) m_conn)->Attribute (name); break;
		default: break;
		}
		return ret;
	}

	
	/**
	 * @brief           Get a text node of an element
	 *
	 *                  @see Configurable::GetText
	 *
	 * @param  path     X-Path
	 * @return          Text
	 */
	inline const char*
	GetText            (const char* path) {
		const char* ret;
		switch (m_ct) {
		case  LOCAL: ret = ( (LocalConnector*) m_conn)->GetText (path); break;
		case REMOTE: ret = ((RemoteConnector*) m_conn)->GetText (path); break;
		default: break;
		}
		return ret;
	}

	
	/**
	 * @brief           Get a text node of an element
	 *
	 *                  @see Configurable::GetText
	 *
	 * @param  path     X-Path
	 * @return          Text
	 */
	inline TiXmlElement*
	GetElement          (const char* path) {
		TiXmlElement* ret;
		switch (m_ct) {
		case  LOCAL: ret = ( (LocalConnector*) m_conn)->GetElement (path); break;
		case REMOTE: ret = ((RemoteConnector*) m_conn)->GetElement (path); break;
		default: break;
		}
		return ret;
	}


private:
	
	ConType m_ct;
	Connection* m_conn; /**< Actual connection */

};
	

}
#endif
