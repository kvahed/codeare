#ifndef __TOOLBOX_HPP__
#define __TOOLBOX_HPP__

#include <vector>
#include <string>

class Toolbox {

public:

	Toolbox() {};
	~Toolbox(){};

	void split (std::vector<std::string>& sv, const std::string& str, const std::string& dlm);

};

#endif
