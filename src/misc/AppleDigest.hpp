#include <CommonCrypto/CommonDigest.h>

inline static char* 
charstar (const std::string& in) {

	size_t len = in.size(); 
	char * buf = new char[len + 1];
	std::copy(in.begin(), in.end(), buf);
	buf[len] = 0;
	return buf;

}

inline static std::string 
sha256 (const std::string& in) {
	
	char* str = charstar (in);
	char buf[65];
    unsigned char hash[CC_SHA256_DIGEST_LENGTH];
    CC_SHA256_CTX sha256;

    CC_SHA256_Init(&sha256);
    CC_SHA256_Update(&sha256, str, strlen(str));
    CC_SHA256_Final(hash, &sha256);
	
    for(size_t i = 0; i < CC_SHA256_DIGEST_LENGTH; i++)
        sprintf(buf + (i * 2), "%02x", hash[i]);
    buf[64] = 0;
	
	delete [] str;
	return std::string(buf);

}
