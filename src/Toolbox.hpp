#ifndef __TOOLBOX_HPP__
#define __TOOLBOX_HPP__

#include <vector>
#include <string>

static std::string bars   = "***************************************************";
static std::string blancs = "                                                   ";

/**
 * @brief  A toolbox for some static stuff
 */
class Toolbox {

public:

	/**
	 * @brief       Singleton destructor
	 */
	~Toolbox();


    /**
     * @brief       Singleton instance
     */
    static Toolbox*  
    Instance        ();


	/**
	 * @brief       String splitter
	 *
	 * @param  sv   Splitted parts
	 * @param  str  String for splitting
	 * @param dlm   Delimiter
	 */
	void 
	Split           (std::vector<std::string>& sv, const std::string& str, const std::string& dlm) const;


	/**
	 * @brief       CPU clock rate.
	 *              Humble abuse of FFTW cycle for timing information.
	 *
	 * @return      clock rate
	 */
	double 
	ClockRate       () const ;	


	
	/**
	 * @brief       Command line progress bar to cout
	 *
	 * @param  pre  Text before bar
	 * @param  post Test after bar
	 * @param  p    Percent
	 */ 
	void 
	ProgressBar    (const std::string& pre, const std::string& post, const short& p) const;

		
private:


	/**
	 * @brief      Hide constructor for singletons
	 */
	Toolbox   () {};

	static Toolbox*    m_instance;           /**< @brief Single instance */

};

#endif
