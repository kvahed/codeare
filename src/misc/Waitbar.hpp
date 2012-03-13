#ifndef __WAITBAR_HPP__
#define __WAITBAR_HPP__

#include <iomanip>
#include <string>

using namespace std;

static string bars   = "***************************************************";
static string blancs = "                                                   ";
static string spear  = ">";

/**
 * @brief       Command line progress bar to cout
 *
 * @param  pre  Text before bar
 * @param  post Test after bar
 * @param  p    Percent
 */ 
static inline void 
ProgressBar    (const string& pre, const string& post, const short& p)  {
	
	assert (p >=  0);
	assert (p <=100);
	
	cout << "\r";
	cout << pre.c_str();
	cout << " | "; 
	cout << bars.substr(0, p/2) << "> " <<  blancs.substr(0, 50-p/2) << "| " << setw(3) << setfill(' ') << p << "% done";
	
}


#endif
