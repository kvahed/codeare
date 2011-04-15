#ifndef __TOOLBOX_HPP__
#define __TOOLBOX_HPP__

#include <vector>
#include <string>

void split (std::vector<std::string>& sv, const std::string& str, const std::string& dlm) {

	assert (dlm.size() > 0);

	size_t  start = 0, end = 0;

	while (end != std::string::npos) {

		end = str.find (dlm, start);

		// If at end, use length=maxLength.  Else use length=end-start.
		sv.push_back(str.substr(start, (end == std::string::npos) ? std::string::npos : end - start));

		// If at end, use start=maxSize.  Else use start=end+delimiter.
		start = ((end > (std::string::npos - dlm.size())) ? std::string::npos : end + dlm.size());
	}
};

#endif
