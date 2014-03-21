#include <openssl/evp.h>


inline static char* 
charstar (const std::string& in) {

	size_t len = in.size(); 
	char * buf = new char[len + 1];
	std::copy(in.begin(), in.end(), buf);
	buf[len] = 0;
	return buf;

}


inline static std::string
digest (const std::string& in, const std::string& type) {

	char* str;
	std::string out;
	char buf[EVP_MAX_MD_SIZE];
	EVP_MD_CTX *mdctx;
	const EVP_MD *md;
	unsigned char md_value[EVP_MAX_MD_SIZE];
	unsigned int md_len;

	OpenSSL_add_all_digests();
	str = charstar(in);
	md = EVP_get_digestbyname(type.c_str());
	mdctx = EVP_MD_CTX_create();
	EVP_DigestInit_ex(mdctx, md, NULL);
	EVP_DigestUpdate(mdctx, str, strlen(str));
	EVP_DigestFinal_ex(mdctx, md_value, &md_len);
	EVP_MD_CTX_destroy(mdctx);

	for(size_t i = 0; i < md_len; i++) 
		sprintf(buf + (i * 2), "%02x", md_value[i]);

	out = std::string (buf); 
	delete[] str;
	return out;

}


inline static std::string 
sha256 (const std::string& in) {
#if defined(_MSC_VER) && (_MSC_VER < 1600)
	return in;
#else
	return digest (in, "sha256");
#endif
}


inline static std::string
sha512 (const std::string& in) {
	return digest (in, "sha512");
}


