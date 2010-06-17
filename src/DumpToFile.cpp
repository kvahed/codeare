#include "DumpToFile.h"

error_code 
DumpToFile::ProcessData () {

	if (m_have_raw) {
		m_raw.dump    ((char*) "raw.bin"   );
		cout << "Dumping raw data" << endl;
	}

	if (m_have_helper)
		m_helper.dump ((char*) "helper.bin");
	if (m_have_pixel)
		m_pixel.dump  ((char*) "pixel.bin" );

	return RRSModule::OK;

}
