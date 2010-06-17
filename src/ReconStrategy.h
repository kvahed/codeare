#ifndef __RECON_STRATEGY_H__
#define __RECON_STRATEGY_H__

#include "Matrix.h"
#include "RRSModule.hh"

#include <cstdlib>
#include <complex>


using namespace std;
using namespace RRSModule;

/**
 * @brief Strategy for reconstruction strategies
 *        Derive hereof to expand the reconstruction toolbox
 *
 */
class ReconStrategy {


public:

	/**
	 * @brief Missing
	 */ 
	ReconStrategy  () {  };
	
	/**
	 * @brief Missing
	 */ 
	virtual
	~ReconStrategy () {};

	/**
	 * @brief Data procession function 
	 */ 
	virtual error_code
	ProcessData     () = 0;
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetRaw           (raw_data* raw)   {
		
		for (int i = 0; i < m_raw.Size(); i++) {
			raw->dabs[i] = abs(m_raw[i]);
			raw->darg[i] = arg(m_raw[i]); 
		}
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetRaw           (const raw_data* raw)   {

		int dim[16];

		m_have_raw = true;

		for (int i = 0; i < 16; i++){
			dim[i] = raw->dims[i];
		}

		m_raw.Reset (dim);

		for (int i = 0; i < m_raw.Size(); i++)
			m_raw[i] =  complex<float> (raw->dabs[i], raw->darg[i]);
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetHelper           (raw_data* helper)   {

		for (int i = 0; i < m_helper.Size(); i++) {
			helper->dabs[i] = abs(m_helper[i]);
			helper->darg[i] = arg(m_helper[i]); 
		}
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetHelper           (const raw_data* helper)   {

		m_have_helper = true;

		int dim[16];

		for (int i = 0; i < 16; i++)
			dim[i] = helper->dims[i];

		m_helper.Reset (dim);

		for (int i = 0; i < m_helper.Size(); i++)
			m_helper[i] =  complex<float> (helper->dabs[i], helper->darg[i]);
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetPixel           (pixel_data* pixel)   {

		for (int i = 0; i < m_pixel.Size(); i++)
			pixel->vals[i] = m_pixel[i];

	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetPixel           (const pixel_data* pixel)   {

		m_have_pixel = true;

		int dim[16];

		for (int i = 0; i < 16; i++)
			dim[i] = pixel->dims[i];

		m_pixel.Reset (dim);

		for (int i = 0; i < m_pixel.Size(); i++) 
			m_pixel[i] =  pixel->vals[i];
		
	};
	
	/**
	 * @brief Get data from recon
	 */
	void 
	GetLabels           (strings labels)   {

		m_labels.length(labels.length());

		for (int i = 0; i < m_labels.length(); i++)
			m_labels[i] = labels[i];
			
	}
	
	/**
	 * @brief Set data for recon
	 */
	void 
	SetLabels          (strings labels)   {

		for (int i = 0; i < m_labels.length(); i++)
			labels[i] = m_labels[i];			

	};
	
	
protected:

	Matrix< complex<float> > m_raw;         /**< raw data matrix   >*/
	Matrix< complex<float> > m_helper;      /**< helper matrix     >*/
	Matrix< short >          m_pixel;       /**< pixel data matrix >*/

	bool                     m_have_raw;
	bool                     m_have_helper;
	bool                     m_have_pixel;

	strings                  m_labels; /**< Labels from the sequence (UID etc)*/

};

#endif /* __RECON_STRATEGY_H__ */
