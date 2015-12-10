/*
 *  codeare Copyright (C) 2007-2015 Kaveh Vahedipour
 *                                  NYU School of Madicine
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
#include <thread>

using namespace codeare::matrix::io;


int main (int argc, char** argv) {

    int error = 1;
	
	if (commandline_opts (argc, argv)) {
        
		// Make local or remote connection
		Connector con (argc, argv, name, debug);
		if (!con.ReadConfig (config_file_uri.c_str())) {
            printf ("*** ERROR: No or corrupt configuration file specified to algorithm! \n");
            return (int)codeare::CONFIG_LOAD_FILE_FAILED;
        }
        
		// Read binary data from file
        TiXmlElement* datain;
        if (infile) { 
            datain = new TiXmlElement("data-in");       // file specified in command line
            datain->SetAttribute("fname", infile);
        } else {                                        
            datain = con.GetElement("/config/data-in"); // file specified in xml
        }

		if (!datain) {
			printf ("*** WARNING: No input data specified to algorithm! \n");
		} else {
			TiXmlElement* datain_entry = datain->FirstChildElement();
			IOContext ic (datain, base_dir, READ);

            if (ic.Strategy()!=SYNGO) {
                while (datain_entry) {
                    
                    const std::string data_name = datain_entry->Value();
                    const std::string data_uri  = datain_entry->Attribute("uri");
                    const std::string data_type = datain_entry->Attribute("dtype");
                    
                    printf ("  Reading %s\n", data_name.c_str());
                    
                    if (data_name.length() == 0) {
                        printf ("*** ERROR: specify a non-empty name for a data set\n");
                        std::cout << "           Entry: " << *datain_entry << std::endl;
                    }
                    
                    if (data_type.length() == 0)
                        printf ("*** ERROR: reading binary data \"%s\" with type %s. Exiting.\n", data_name.c_str(), data_type.c_str());
                    
                    // TODO: check first if entry exists and has right format
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
                    } else  {
                        printf ("*** ERROR: Couldn't load a data set specified in\n");
                        std::cout << "           Entry: " << *datain_entry << std::endl;
                        return codeare::WRONG_OR_NO_DATASET;
                    }
                    
                    datain_entry = datain_entry->NextSiblingElement();
                }
            } else {
            	ic.Read();
            }
        }

		TiXmlElement* chain = con.GetElement("/config/chain");
		if (!chain) {
			printf ("*** ERROR: Chain of modules must be specified. \n");
			return codeare::CONFIG_MISSING_CHAIN;
		}
	    TiXmlElement* module = chain->FirstChildElement();
	    if (!module) {
			printf ("*** ERROR: Chain of modules must have at least one element.\n");
			return codeare::CONFIG_EMPTY_CHAIN;
		}

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
	    	codeare::error_code ce = con.Init (module_name.c_str(), config.c_str());
	    	if (ce != codeare::OK) {
	    		printf ("  *** ERROR: Intialising failed ... bailing out!\n");
	    		return (int) ce;
	    	}
	    	++nmodules;
	    	module = module->NextSiblingElement();

	    }

	    if (nmodules > 0) {
	    	con.Prepare();
	    	con.Process();
	    } else {
	    	printf ("Warning! No modules were found in the configuration file. Exiting\n");
	    }

		TiXmlElement* dataout = con.GetElement("/config/data-out");
		if (!dataout)
			printf ("*** WARNING: No output data expected from algorithm? \n");

		else {

			TiXmlElement* dataout_entry = dataout->FirstChildElement();
			IOContext out = IOContext (dataout, base_dir, WRITE);

			while (dataout_entry) {
				const std::string data_name = dataout_entry->Value();
				const std::string data_type = dataout_entry->Attribute("dtype");
				codeare::error_code ec;

				if (!(data_name.length() && data_type.length()))
					printf("Error writing binary data \"%s\" with type %s. Exiting.\n", data_name.c_str(), data_type.c_str());

				if        (TypeTraits<float>::Abbrev().compare(data_type) == 0)  {
					Matrix<float> M;
					if ((ec = con.GetMatrix(data_name, M)) == codeare::OK)
						out.Write(M, dataout_entry);
				} else if (TypeTraits<double>::Abbrev().compare(data_type) == 0) {
					Matrix<double> M;
					if ((ec = con.GetMatrix(data_name, M)) == codeare::OK)
						out.Write(M, dataout_entry);
				} else if (TypeTraits<cxfl>::Abbrev().compare(data_type) == 0)   {
					Matrix<cxfl> M;
					if ((ec = con.GetMatrix(data_name, M)) == codeare::OK)
						out.Write(M, dataout_entry);
				} else if (TypeTraits<cxdb>::Abbrev().compare(data_type) == 0)   {
					Matrix<cxdb> M;
					if ((ec = con.GetMatrix(data_name, M)) == codeare::OK)
						out.Write(M, dataout_entry);
				}

				if (ec == codeare::NO_MATRIX_IN_WORKSPACE_BY_NAME) {
					printf ("*** ERROR: No dataset by name \"%s\" found ... bailing out!\n", data_name.c_str());
					return (int) ec;
				} else if (ec == codeare::WRONG_MATRIX_TYPE) {
					printf ("*** ERROR: Dataset by name \"%s\" found. Wrong type requested though ... bailing out!\n", data_name.c_str());
					return (int) ec;
				}

				dataout_entry = dataout_entry->NextSiblingElement();
			}
		}
		
	    con.Finalise();

        error = 0;

	} 
		
    return error;	

}



