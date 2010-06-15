#ifndef __RECONCONTEXT_H__
#define __RECONCONTEXT_H__

#include <vector>
#include "ReconStrategy.h"

using namespace std;

/**
 * @brief Context of a reconstruction method
 */
class ReconContext {
  

  
public:

	/**
	 * @brief Construct with setting up available strategies
	 */
	ReconContext ();

	/**
	 * @brief Construct with a strategy
	 */
	ReconContext (ReconStrategy* strategy) {};

	/**
	 * @brief Alter strategy
	 */
	inline void 
	Strategy     (ReconStrategy* strategy) {
		m_strategy = strategy;
	}
	
	/**
	 * @brief Alter strategy
	 */
	inline void 
	Strategy     (method m) {
		m_strategy = m_strategies.at((int) m);
	}
	
	/**
	 * @brief get active startegy
	 */
	inline ReconStrategy*
	Strategy     () {
		return m_strategy;
	}
	
	/**
	 * @brief Process data with given strategy
	 */
	inline RRSModule::error_code
	ProcessData (method m) {
		return m_strategy->ProcessData();
	}
	

private:

	ReconStrategy*            m_strategy;   /**< Active strategy      */
	vector< ReconStrategy* >  m_strategies;  /**< Available strategies */
	
};

#endif 
