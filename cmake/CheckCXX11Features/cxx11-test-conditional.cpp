#include <type_traits>
#include <string>

int main() {
	typedef typename std::conditional<false, const std::string, std::string>::type StringType;
    StringType s;
	return 0;
}
