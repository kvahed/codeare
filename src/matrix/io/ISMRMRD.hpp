/*
 * ISMRMRD.hpp
 *
 *  Created on: Jan 8, 2013
 *      Author: kvahed
 */

#ifndef __ISMRMRD_HPP__
#define __ISMRMRD_HPP__

#include "ismrmrd_hdf5.h"
#include "ismrmrd.hxx"

namespace matrix {

	namespace io {

		namespace ismrm {

			class IRDFile {

			public:

				/**
				 * @brief  Construct with file name
				 *
				 * @param  fname    File name
				 * @param  verbose  Verbose (default: false)?
				 */
				IRDFile (const std::string& fname, const bool verbose = false) {
					m_props.schema_location ("http://www.ismrm.org/ISMRMRD", m_scheme);
					m_fname = fname;
				}

				/**
				 * @brief  Clean up
				 */
				~IRDFile () {}

				bool IORead (const std::string& dname, const bool verbose = false) {

					ISMRMRD::IsmrmrdDataset ds (m_fname.c_str(),dname.c_str());

					size_t i = 0, na = ds.getNumberOfAcquisitions();

					// Loop over acquisitions
					for (; i < na; i++) {

						boost::shared_ptr<ISMRMRD::Acquisition> acq = ds.readAcquisition(i);
						//size_t offset = acq->head_.idx.kspace_encode_step_1*buffer.dimensions_[0];
						//memcpy(&buffer.data_[offset],acq->data_,sizeof(float)*2*buffer.dimensions_[0]);

					}

					return true;

				}

				xml_schema::properties m_props; /**< @brief Properties */

				std::string m_scheme;           /**< @brief Scheme */
				std::string m_fname;

			};

		}

	};

};


#endif /* __ISMRMRD_HPP__ */
