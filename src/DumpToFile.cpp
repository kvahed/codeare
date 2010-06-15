#include "DumpToFile.h"

error_code 
DumpToFile::ProcessData () {

	cout << m_labels[0] << endl;

	//filename = "out.bin";

	//for (int i=0; i++)
	//ofstream fout(fname.c_str() , ios::binary);
	//fout.write((char*)(&(mag)), sizeof(mag));
	//fout.write((char*)(&(pha)), sizeof(pha));

	return RRSModule::OK;

}
