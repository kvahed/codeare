#if defined(__APPLE__)
#  include "AppleDigest.hpp"
#else
#  include "Digest.hpp"
#endif

template <class T> bool 
sha256test (Connector<T>* rc) {

	std::string input;
	std::string coded;
	std::string hw = "Hello World";

	std::string shao;

	shao = sha256 (hw);
	printf ("\n'Hello World' - %s\n\n", shao.c_str());
	
	return true;

}
