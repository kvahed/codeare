/*
 * SyngoMRProtocol.h
 *
 *      Author: kvahed
 *      Project: codeare
 */

#ifndef _SYNGO_MR_PROTOCOL_H_
#define _SYNGO_MR_PROTOCOL_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/regex.hpp>

#include <string>
#include <fstream>
#include <iostream>

#if (Boost_MINOR_VERSION>55)
typedef boost::property_tree::ptree::key_type settings_key_t;
#else
typedef char settings_key_t;
#endif


/**
 * @brief Type conversion from string
 */
template <typename T> struct ConvertTo;

template<> struct ConvertTo<short> {
	static short str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<int> {
	static int str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<long> {
	static long str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<unsigned short> {
	static unsigned short str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<unsigned int> {
	static unsigned int str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<unsigned long> {
	static unsigned long str (const std::string& in) { return std::stoi(in); }
};
template<> struct ConvertTo<float> {
	static float str (const std::string& in) { return std::stof(in); }
};
template<> struct ConvertTo<double> {
	static double str (const std::string& in) { return std::stof(in); }
};
template<> struct ConvertTo<const char*> {
	static const char* str (const std::string& in) { return in.c_str(); } 
};
template<> struct ConvertTo<std::string> {
	static std::string str (const std::string& in) { return std::string(in); } 
};

/**
 * @brief Search result
 **/
typedef struct search_result {
	long _pos;
	std::string _str;

	/**
	 * @brief Construct with position and found string
	 * @param pos position
	 * @param str found string
	 */
	search_result (long pos = LONG_MAX, const std::string& str = "") :
		_pos(pos), _str(str) {};

	/**
	 * @brief Get position
	 * @return position in search text
	 */
	inline long position () const { return _pos;}

	/**
	 * @brief Get found string
	 * @return found string
	 */
	inline const std::string& str () const { return _str; }

	/**
	 * @brief Get length of found string
	 * @return length of found string
	 */
	inline long length () const { return _str.length(); }

	/**
	 * @brief Empty?
	 * @return empty search?
	 */
	inline bool empty() const { return _str.empty(); }

	/**
	 * @brief Found?
	 * @return found anything?
	 */
	inline bool found() const { return !_str.empty(); }

	/**
	 * @brief Dump found string to output stream
	 * @param output stream to dump to
	 * @param the search result
	 * @return the output stream to dump to
	 */
	friend std::ostream& operator<< (std::ostream& os,
			const search_result& res) {	os << res.str(); return os;	}

} search_result;

/**
 * @brief Protocol representation Syngo MR
 */
class SyngoMRProtocol {

public:

	typedef boost::property_tree::ptree stack_item;

	enum read_exception {
		CANNOT_OPEN = 100,
		ZERO_HEADER_LENGTH,
		TOO_LARGE_HEADER
	};

	enum parse_exception {
		FAILED_TO_FIND_XPROTOCOL = 200,
		HIT_EOF_BEFORE_FINISHING,
		PROTOCOL_HIERARCHY_ERROR
	};

	enum get_exception {
		INVALID_PATH = 300,
		CONVERSION_FAILED
	};

	enum out_type {
		PLAIN, XML, JSON, INI
	};

	enum major_version {VAB, VD};

	/**
	 * @default Constructor
	 */
	SyngoMRProtocol();

	/**
	 * @brief Construct with file name
	 */
	SyngoMRProtocol (const std::string& filename, int verbosity = 0, bool remove_empty_tags = true);

	/**
	 * @brief Copy
	 */
	SyngoMRProtocol (const SyngoMRProtocol& protocol);

	/**
	 * @brief Clean up
	 */
	virtual ~SyngoMRProtocol();

	/**
	 * @brief Property object
	 */
	virtual const boost::property_tree::ptree& Properties () const;

	/**
	 * @brief Get file name
	 */
	virtual const std::string& FileName () const;

	/**
	 * @brief Raw output
	 */
	virtual const std::string& Raw () const;

	/**
	 * @brief XML output
	 */
	virtual std::ostream& ToXML (std::ostream& os) const;
	virtual void ToXML (const std::string& filename) const;

	/**
	 * @brief Get a value from property tree
	 */
	const std::string GetStr (const std::string& path) const;

	template<typename T> T Get (const std::string& path) const {
		T ret = 0;
		std::string str_rep;
		try {
			str_rep = GetStr(path);
            ret = ConvertTo<T>::str(str_rep);
		} catch (const boost::property_tree::ptree_bad_path&) {
			printf ("  ERROR: Invalid path %s !", path.c_str());
			throw INVALID_PATH;
		} catch (const std::logic_error& e ) {
		} catch (const std::exception&) {
			printf ("  ERROR: Conversion to %s failed for %s !",
					typeid(T).name(), str_rep.c_str());
			throw CONVERSION_FAILED;
		}
		return ret;
	}

	friend std::ostream& operator<< (std::ostream &out, const SyngoMRProtocol &smrp);

	/**
	 * @brief Get the offset of the last measurement in the raw file
	 */
	uint64_t GetLastMeaseOffset() const;

private:

	/**
	 * @brief Reads file into RAM buffer for parsing
	 * @param level of verbosity
	 * @return happiness
	 */
	virtual bool ReadFile (int verbosity = 0);

	/**
	 * @brief Parse RAM buffer for content
	 * @param level of verbosity
	 * @return happiness
	 */
	virtual bool Parse (int verbosity = 0);

	virtual void TrimData ();
	//virtual void GotoXProtocol (const std::string& section = std::string("Config"));
	virtual std::string GetNextLine (size_t pos = 0);

	void HandleOpenBracket (std::vector<int>& non_stack, long pos);
	void HandleCloseBracket (std::vector<stack_item*>& stack,
			std::vector<int>& non_stack, long pos);
	void HandleNode (std::vector<stack_item*>& stack, std::vector<int>& non_stack);
    void HandleAscConvEntry (stack_item& root, const std::string& ascconv_str);

	std::string _raw;  // measurement header
	std::string _meas_file_name; // measurement file name
	uint32_t _header_len;
	uint32_t _data_len;
	major_version _mver;
	std::string::const_iterator _xpos;
	boost::property_tree::ptree _props;
	std::string::const_iterator _start, _cur, _end;

	// A flag indicating if the empty tags shall be removed when doing TrimData()
	uint64_t _last_meas_offset;
	bool _remove_empty_tags;
};

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>

inline static const std::string RegexReplace (const std::string& str,
		const std::string& what, const std::string& by = "") {
	std::string ret = boost::regex_replace(str, boost::regex(what), by);
	return ret;
}

inline SyngoMRProtocol::SyngoMRProtocol (const std::string& file_name, int verbosity, bool remove_empty_tags) :
	_meas_file_name (file_name), _last_meas_offset(0), _remove_empty_tags(remove_empty_tags) {

	if (verbosity > 0)
		printf ("Syngo MR protocol \n  File name: %s\n", _meas_file_name.c_str());

	try {
		ReadFile (verbosity);
		try {
			Parse(verbosity);
		} catch (const SyngoMRProtocol::parse_exception& status) {
			switch (status) {
			case FAILED_TO_FIND_XPROTOCOL:
				printf ("  ERROR: Could not find the XProtocol.\n");
				break;
			case HIT_EOF_BEFORE_FINISHING:
				printf ("  ERROR: Hit EOF before finishing parsing.\n");
				break;
			case PROTOCOL_HIERARCHY_ERROR:
				printf ("  ERROR: Protocol hierarchy error.\n");
				break;
			default:
				assert(false); // We should have never get here
				break;
			}
		}
	} catch (const SyngoMRProtocol::read_exception& status) {
		switch (status) {
		case CANNOT_OPEN:
			printf ("  ERROR: Filed to open file %s.\n", _meas_file_name.c_str());
			break;
		case TOO_LARGE_HEADER:
			printf ("  ERROR: Header way too long %d.\n", _header_len);
			break;
		case ZERO_HEADER_LENGTH:
			printf ("  ERROR: Zero header length.\n");
			break;
		default:
			assert(false); // We should have never get here
			break;
		}
		printf ("  EXITING!\n");
	}

}

inline SyngoMRProtocol::SyngoMRProtocol() : _data_len(0), _last_meas_offset(0),
		_mver(VAB), _remove_empty_tags(true), _header_len(0) {}

inline SyngoMRProtocol::~SyngoMRProtocol() {
	_raw.clear();
}

inline bool SyngoMRProtocol::ReadFile (int verbosity) {

	std::ifstream file; // file handle
	using namespace codeare::matrix::io;
	// Open
	if (verbosity > 0)
		printf ("  Opening ...");
	file.open (_meas_file_name.c_str());
	if (!file.is_open())
		throw CANNOT_OPEN;
	if (verbosity > 0)
		printf (" done.\n");

	// VA/VB or VD?
	uint32_t x[2];
	file.read ((char*)x, 2*sizeof(uint32_t));
	if (x[0] == 0 && x[1] <= 64)
		_mver = VD;
	else
		_mver = VAB;
	file.seekg (0);

	if (_mver == VD) {
        uint32_t id, ndset;
        std::vector<VD::EntryHeader> veh;
        file.read ((char*)&id, sizeof(uint32_t));      // ID
        file.read ((char*)&ndset, sizeof(uint32_t));   // # data sets
        if (verbosity > 0)
        	printf ("    VD Raid file (id:%d) contains %d data set(s):\n", id, ndset);
        veh.resize(ndset);
        for (size_t i = 0; i < ndset; ++i) {
            file.read ((char*)&veh[i], VD::ENTRY_HEADER_LEN);
            if (verbosity > 0)
            	printf ("       MID(%d) FID(%d) %s \n", veh[i].MeasID, veh[i].FieldID, veh[i].ProtocolName);
        }

		_last_meas_offset = veh.back().MeasOffset;
        file.seekg(_last_meas_offset);                // Go to last measurement
	}

	// Find header length
	file.read((char*)&_header_len, sizeof(uint32_t));
	if (_header_len <= 0)
		throw ZERO_HEADER_LENGTH;
	else if (_header_len > 1000000)
		throw TOO_LARGE_HEADER;

	// Read header
	if (verbosity > 0)
		printf ("    Header holds %d bytes.\n  Reading ...", _header_len);
	file.seekg(0);
	_raw.resize(_header_len);
	file.read(&_raw[0], _header_len);
	_start = _raw.begin();
	_cur = _start;
	_end = _raw.end();
	if (verbosity > 0)
		printf (" done.\n");

	// Close file
	if (verbosity > 0)
		printf ("  Closing ...");
	file.close ();
	if (verbosity > 0)
		printf (" done.\n");

	return true;

}

inline std::string SyngoMRProtocol::GetNextLine (size_t pos) {
	return _raw.substr(pos,_raw.find("<>", pos));
}

inline void SyngoMRProtocol::TrimData () { // Thanks, Siemens!
	_raw   = RegexReplace (_raw, "<Default>\\s+<", "<");
	_raw   = RegexReplace (_raw, "<ParamDouble.\"\"> \\{ <Precision> 16 \\}"
			"|<Label>\\s*[\\s|a-z|A-Z|0-9|\\.|\"]+"
		 	"|<Precision> [0-9]*|<M(in|ax)Size> [0-9]*|<Tooltip> \"*\""
			"|<Default>\\s*[0-9|a-z|A-Z|\\.|\"|\\s|\\-]+\\s*");
	if (_remove_empty_tags) {
		_raw = RegexReplace(_raw, "\\{\\s*\\}");
	}
	_raw   = RegexReplace (_raw, "<Limit>\\s*[0-9|a-z|A-Z|\\.|\"|\\s]+\\s*");
	if (_remove_empty_tags) {
		_raw = RegexReplace(_raw, "\\{\\s*\\}"
				"|<[a-z|A-Z|\"|\\.]+>\\s*\\{\\s*\\s*\\}");
	}
	_raw   = RegexReplace (_raw, "<LimitRange>\\s*[\\s|a-z|A-Z|0-9|\\.|\"]+");
	if (_remove_empty_tags) {
		_raw = RegexReplace(_raw, "<[a-z|A-Z|\"|\\.|_]+>\\s*\\{\\s*\\s*\\}");
		_raw = RegexReplace(_raw, "<[a-z|A-Z|\"|\\.|_]+>\\s*\\{\\s*\\s*\\}");
		_raw = RegexReplace(_raw, "<[a-z|A-Z|\"|\\.|_]+>\\s*\\{\\s*\\s*\\}");
		_raw = RegexReplace(_raw, "<[a-z|A-Z|\"|\\.|_]+>\\s*\\{\\s*\\s*\\}");
	}
	_raw   = RegexReplace (_raw, "\\.\"1", "\\.\"One$1");
	_raw   = RegexReplace (_raw, "\\.\"2", "\\.\"Two$1");
	if (_remove_empty_tags) {
		_raw   = RegexReplace (_raw, "(<[a-z|A-Z]*\\.\"[a-z|A-Z]+\">\\s*)+\\}", "\\}");
		_raw   = RegexReplace (_raw, "<[a-z|A-Z]*\\.\"[a-z|A-Z]*\">\\s*\\{\\s*\\}");
	}
	_raw   = RegexReplace (_raw, "[<[a-z|A-Z]*\\.\"[a-z|A-Z]*\">\\s+]+");
	_raw   = RegexReplace (_raw, "\\s*\\}\\s*\\{\\s*", ", ");
	_raw   = RegexReplace (_raw, "\\{\\s*, ", "\\{ ");
	_raw   = RegexReplace (_raw, "(\\-[1|2]\\s*)+");
	_raw   = RegexReplace (_raw, ">\\s+\\{", "> \\{");
	_raw   = RegexReplace (_raw, "(\\s*\\n+)+", "\\n");
	_raw   = RegexReplace (_raw, "(<[a-z|A-Z]+\\.*\"*[a-z|A-Z]*\"*>\\s*)+<", "<");
	_raw   = RegexReplace (_raw, "<Class> \"[a-z|A-Z|\\@|_]*\"");
	_raw   = RegexReplace (_raw, "\\{\\s*\\{\\s*([a-z|A-Z|0-9|\"|\\.]+)\\s*\\}\\s*\\}", "\\{$1\\}");
	_raw   = RegexReplace (_raw, "\\{\\s*\\{\\s*([a-z|A-Z|0-9|\"|\\.]+)\\s*\\}\\s*\\}", "\\{$1\\}");

	_raw   = RegexReplace (_raw, "([A-Z][a-zA-Z]+).{3,10}<XProtocol>", "<$1>");
	_start = _raw.begin();
	_cur   = _start;
	_end   = _start + _raw.size();
}

inline void SyngoMRProtocol::HandleOpenBracket (std::vector<int>& non_stack, long pos) {
	non_stack.push_back(0);
	_cur += pos+1;
}

inline void SyngoMRProtocol::HandleCloseBracket (std::vector<stack_item*>& stack,
		std::vector<int>& non_stack, long pos) {
	if (non_stack.size()>0)
		non_stack.pop_back();
	else if (stack.size()>1)
		stack.pop_back();
	_cur += pos+1;
}

inline static search_result RegexSearch (std::string::const_iterator from,
		std::string::const_iterator to,	const std::string& what) {
	boost::match_results<std::string::const_iterator> mr;
	if (boost::regex_search (from, to, mr, boost::regex(what)))
		return search_result(mr.position(), mr.str());
	return search_result();
}

inline void SyngoMRProtocol::HandleNode (std::vector<stack_item*>& stack, std::vector<int>& non_stack) {
	static std::string node_regexp = "<([a-z|A-Z]*)\\.*\"*([a-z|A-Z|0-9|_]*)\"*>\\s*\\{[a-z|A-Z|0-9|\"|;|\\\\|\\%|\\.|_|/|,|\\-|\\&|\\{|\\}|\\s]*",
			leaf_regexp = "<([a-z|A-Z]*)\\.*\"*([a-z|A-Z|0-9|_]*)\"*>\\s*\\{([a-z|A-Z|0-9|\"|;|\\\\|\\%|\\.|_|/|,|\\-|\\&|\\s]+)\\}";
	search_result res = RegexSearch (_cur, _end, node_regexp);

	_cur += res.position();
	std::string node_name, node_type, node_val;
	if (res.found()) {
		search_result leaf = RegexSearch (res.str().begin(), res.str().end(), leaf_regexp);
		if (leaf.found()) { 	// leaf
			_cur += leaf.length();
			node_type = RegexReplace (leaf.str(), leaf_regexp, "$1");
			node_name = RegexReplace (leaf.str(), leaf_regexp, "$2");
			node_val  = RegexReplace (leaf.str(), leaf_regexp, "$3");
			if (node_name == "")
				node_name = node_type;
			stack.back()->put (node_name, node_val);
		} else {
			_cur += res.length();
			node_type = RegexReplace (res.str(), node_regexp, "$1");
			node_name = RegexReplace (res.str(), node_regexp, "$2");
			if (node_name == "")
				node_name = node_type;
			search_result node = RegexSearch (res.str().begin(), res.str().end(), node_regexp);
			if (node_name == "Dicom" || node_name == "Meas" ||
					node_name == "Phoenix"|| node_name == "Spice") {
				non_stack.clear();
				boost::property_tree::ptree* root = stack[0];
				stack.clear();
				stack.push_back(root);
			}
			stack.push_back(&stack.back()->put (node_name, ""));
		}
	}

}

inline void SyngoMRProtocol::HandleAscConvEntry (boost::property_tree::ptree& root, const std::string& ascconv_str) {
    boost::property_tree::ptree& ascconv = root.put("ASCCONV", "");
    std::istringstream istr (ascconv_str);
    std::string line;
    for (std::string line; std::getline(istr, line); ) {
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
        auto sind = line.find("=");
        std::string key = line.substr(0, sind), val = line.substr(sind+1, line.size()-sind+1);
        try {
            float f = std::stof(val);
            if (std::round(f)==f)
                ascconv.put (key, std::lround(f));
            else
                ascconv.put (key, f);
        } catch (const std::invalid_argument&) {
            ascconv.put (key, val);
        }
    }
}


inline bool SyngoMRProtocol::Parse (int verbosity) {

	// Keep track of tree structure
	std::vector<stack_item*> stack;
	std::vector<int> non_stack;
	boost::property_tree::ptree& root = _props.put ("XProtocol", "");
    stack.push_back(&root);

    // ASCCONV
    size_t abeg = _raw.find("### ASCCONV BEGIN ###")+22;
    HandleAscConvEntry(root, _raw.substr(abeg, _raw.find("### ASCCONV END ###",abeg)-1-abeg));

    // XProtocol
	if (verbosity > 0)
		printf ("  Parsing ...");

	TrimData(); // Get rid of crap

	// Actually parse
 	while (_cur < _end) {

		long dpos = RegexSearch( _cur, _end, "<").position(),
		     opos = RegexSearch( _cur, _end, "\\{").position(),
			 cpos = RegexSearch( _cur, _end, "\\}").position();
		long next = std::min(dpos, std::min(opos, cpos));

		if (opos == next)             // Keep track of open curly
			HandleOpenBracket (non_stack, opos);
		else if (cpos == next)      // Close one curly
			HandleCloseBracket (stack, non_stack, cpos);
		else if (dpos == next)
			HandleNode (stack, non_stack);
		else if (next == LONG_MAX)
			break;
		else
			break;

	}


	if (verbosity > 0)
		printf (" done.\n");

	return true;
}

inline const std::string SyngoMRProtocol::GetStr (const std::string& path) const {
	return RegexReplace(_props.get<std::string>(path),
			"\\s*\"{1}\\s*([a-z|A-Z|0-9|\"|\\.|_|/|,|\\-|\\{|\\}|\\s]*)\\s*\"\{1}\\s*", "$1");
}


inline const boost::property_tree::ptree& SyngoMRProtocol::Properties () const {
	return _props;
}

inline const std::string& SyngoMRProtocol::FileName () const {
	return _meas_file_name;
}

inline const std::string& SyngoMRProtocol::Raw () const {
	return _raw;
}

inline std::ostream& SyngoMRProtocol::ToXML (std::ostream& os) const {
	boost::property_tree::xml_writer_settings<settings_key_t> settings(' ', 2u);
	write_xml(os, _props, settings);
	return os;
};

inline void SyngoMRProtocol::ToXML (const std::string& filename) const {
	boost::property_tree::xml_writer_settings<settings_key_t> settings(' ', 2u);
	write_xml(filename, _props, std::locale(), settings);
};

inline std::ostream& operator<< (std::ostream& os, const SyngoMRProtocol& smrp) {
	return os << smrp.Raw();
}

inline uint64_t SyngoMRProtocol::GetLastMeaseOffset() const
{
	return _last_meas_offset;
}


#endif /* _SYNGO_MR_PROTOCOL_H_ */
