#include "DumpToFile.h"

error_code 
DumpToFile::ProcessData () {

	printf ("Dumping ... ");

	if (m_have_raw) {
		printf ("raw data.\n");
		m_raw.dump    ((char*) "raw.bin"   );
	}

	if (m_have_helper) {
		printf ("helper data.\n");
		m_helper.dump ((char*) "helper.bin");
	}

	if (m_have_pixel) {
		printf ("pixel data.\n");
		m_pixel.dump  ((char*) "pixel.bin" );
	}

	if (m_labels.length() > 0) {
		printf ("pixel data.\n");

	    ofstream of;

	    of.open ("labels.txt");
		for (int i = 0; i < m_labels.length(); i++)
			of << m_labels[i] << "\n";
	    of.close();

	}

	return RRSModule::OK;

}
