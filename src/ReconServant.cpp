/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Jülich, Germany
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

#include "ReconServant.h"
#include "ReconContext.h"


/**************************************************************************************************/
ReconServant::ReconServant  ()               {

	m_config = new char;

}


/**************************************************************************************************/
ReconServant::~ReconServant ()               {

	delete m_config;

}


/**************************************************************************************************/
error_code
ReconServant::Process  (const char* name)       {

	error_code e = OK;

	cout << "Setting incoming data ... " << endl;

	ReconContext* context = new ReconContext(name);

	context->Strategy()->SetRaw(&m_raw);
	context->Strategy()->SetHelper(&m_helper);
	context->Strategy()->SetPixel(&m_pixel);

	context->Strategy()->SetConfig(m_config);
	
	cout << "... done. Will invoke data procession ... " << endl;
	
	e = context->Strategy()->Process();
	
	cout << "... done. Getting processed data..." << endl;
	
	context->Strategy()->GetRaw(&m_raw);
	context->Strategy()->GetHelper(&m_helper);
	context->Strategy()->GetPixel(&m_pixel);
	context->Strategy()->GetConfig(m_config);
	
	cout << "... done. Deleting context ..." << endl;

	delete context;

	cout << "... done. Will handle control back to client." << endl;
	
	return e;
	
}


/**************************************************************************************************/
void
ReconServant::raw          (const raw_data& d)   {
	m_raw = d;
}


/**************************************************************************************************/
raw_data*
ReconServant::raw          ()                    {
	return new raw_data (m_raw);
}


/**************************************************************************************************/
void
ReconServant::helper       (const raw_data& d)   {
	m_helper = d;
}


/**************************************************************************************************/
raw_data*
ReconServant::helper       ()                    {
	return new raw_data (m_helper);
}


/**************************************************************************************************/
void
ReconServant::pixel        (const pixel_data& d) {
	m_pixel = d;
}


/**************************************************************************************************/
pixel_data*
ReconServant::pixel        ()                    {
	return new pixel_data (m_pixel);
}


/**************************************************************************************************/
void 
ReconServant::config       (const char* d)    {
	stringstream tmp;
	tmp << d;
	cout << tmp.str().length() << endl;
	m_config = new char[tmp.str().length()*sizeof(char*)];
	strcpy (m_config, tmp.str().c_str());
}


/**************************************************************************************************/
char* 
ReconServant::config       ()                    {
	return CORBA::string_dup(m_config);
}

