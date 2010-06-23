#include "DumpToFile.h"

error_code 
DumpToFile::ProcessData () {

	if (m_have_raw)
		m_raw.dump    ((char*) "raw.bin"   );
	if (m_have_helper)
		m_helper.dump ((char*) "helper.bin");
	if (m_have_pixel)
		m_pixel.dump  ((char*) "pixel.bin" );

	if (m_labels.length() > 0) {
	    ofstream of;
	    of.open ("labels.txt");
		for (int i = 0; i < m_labels.length(); i++)
			of << m_labels[i] << "\n";
	    of.close();
	}

	return RRSModule::OK;

}
