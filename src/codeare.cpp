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

#include "codeare.hpp"
#include "IOContext.hpp"

using namespace codeare::matrix::io;

int main (int argc, char** argv) {

    int error = 1;
	
	if (commandline_opts (argc, argv)) {

		// Make local or remote connection
		Connector con (argc, argv, name, verbose);
		con.ReadConfig (config_file_uri.c_str());

		TiXmlElement* datain = con.GetElement("/config/data-in");
	    TiXmlElement* datain_entry = datain->FirstChildElement("item");
	    IOContext ic (datain, base_dir, READ);

	    while (datain_entry) {
	    	const std::string data_name = datain_entry->Attribute("uri");
	    	const std::string data_type = datain_entry->Attribute("dtype");
	    	if (!(data_name.length() && data_type.length()))
	    		printf("Error reading binary data %s with type %s. Exiting.\n", data_name.c_str(), data_type.c_str());
	    	if        (TypeTraits<float>::Abbrev().compare(data_type) == 0)  {
	    		Matrix<float> M = ic.Read<float>(datain_entry);
	    		con.SetMatrix(data_name, M);
	    	} else if (TypeTraits<double>::Abbrev().compare(data_type) == 0) {
	    		Matrix<double> M = ic.Read<double>(datain_entry);
	    		con.SetMatrix(data_name, M);
	    	} else if (TypeTraits<cxfl>::Abbrev().compare(data_type) == 0)   {
	    		Matrix<cxfl> M = ic.Read<cxfl>(datain_entry);
	    		con.SetMatrix(data_name, M);
	    	} else if (TypeTraits<cxdb>::Abbrev().compare(data_type) == 0)   {
	    		Matrix<cxdb> M = ic.Read<cxdb>(datain_entry);
	    		con.SetMatrix(data_name, M);
	    	}
	    	datain_entry = datain_entry->NextSiblingElement();
	    }

	    TiXmlElement* chain = con.GetElement("/config/chain");
	    TiXmlElement* module = chain->FirstChildElement();
	    size_t nmodules = 0;

	    while (module) {

	    	std::string config;
	    	config << *module;
	    	std::string module_name = module->Value();
	    	if (module_name.length() == 0) {
	    		printf ("  *** ERROR: Module has no name \"%s\" ... bailing out\n", module_name.c_str());
	    		nmodules = 0;
	    		break;
	    	}

	    	printf ("Initialising %s ...\n", module_name.c_str());
	    	if (con.Init (module_name.c_str(), config.c_str()) != codeare::OK) {
	    		printf ("  *** ERROR: Intialising failed ... bailing out!");
	    		return false;
	    	}
	    	module = module->NextSiblingElement();
	    	nmodules++;
	    }


	    if (nmodules > 0) {
	    	cout << "OOps" << endl;
	    	con.Prepare();
	    	cout << "OOps" << endl;
	    	con.Process();
	    	cout << "OOps" << endl;

	    } else {
	    	printf ("Warning! No modules were found in the configuration file. Exiting\n");
	    }

	    con.Finalise();

        error = 0;

	} 
		
    return error;	

}



