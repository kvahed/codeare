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
#include "IOFile.hpp"

namespace codeare {

	namespace matrix {

		namespace io {

			class IRDFile : public IOFile {

			public:

				/**
				 * @brief  Construct with file name
				 *
				 * @param  fname    File name
				 * @param  verbose  Verbose (default: false)?
				 */
				IRDFile (const std::string& fname,
						const IOMode& mode = READ,
						Params params = Params(),
						const bool verbose = false) :
						IOFile(fname, mode, params, verbose) {

					const char* sl = m_params.Get<std::string>("scheme").c_str();
					m_props.schema_location ("http://www.ismrm.org/ISMRMRD", sl);

				}

				/**
				 * @brief  Clean up
				 */
				~IRDFile () {}

				template<class T> Matrix<T> Read (const std::string& dname) const {

					Matrix<T> M;

					ISMRMRD::IsmrmrdDataset ds (m_fname.c_str(),dname.c_str());

					size_t i = 0, na = ds.getNumberOfAcquisitions();

					boost::shared_ptr<std::string> xml = ds.readHeader();
					std::istringstream str_stream(*xml, std::stringstream::in);

					//Let's print some information from the header
					boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg;

					try {
						cfg = boost::shared_ptr<ISMRMRD::ismrmrdHeader>(ISMRMRD::ismrmrdHeader_ (str_stream,0,m_props));
					}  catch (const xml_schema::exception& e) {
						std::cout << "Failed to parse XML Parameters: " << e.what() << std::endl;
					}

					ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
					if (e_seq.size() != 1) {
						std::cout << "Number of encoding spaces: " << e_seq.size() << std::endl;
						std::cout << "This simple reconstruction application only supports one encoding space" << std::endl;
					}

					ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
					ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
					ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

					std::cout << "Encoding Matrix Size        : [" << e_space.matrixSize().x() << ", " << e_space.matrixSize().y() << ", " << e_space.matrixSize().z() << "]" << std::endl;
					std::cout << "Reconstruction Matrix Size  : [" << r_space.matrixSize().x() << ", " << r_space.matrixSize().y() << ", " << r_space.matrixSize().z() << "]" << std::endl;
					std::cout << "Number of acquisitions      : " << ds.getNumberOfAcquisitions() << std::endl;

					if (e_space.matrixSize().z() != 1) {
						std::cout << "This simple reconstruction application only supports 2D encoding spaces" << std::endl;
					}


					// Loop over acquisitions
					for (; i < na; i++) {

						boost::shared_ptr<ISMRMRD::Acquisition> acq = ds.readAcquisition(i);

						//acq->head_.idx.kspace_encode_step_1
						//size_t offset = acq->head_.idx.kspace_encode_step_1*buffer.dimensions_[0];
						//memcpy(&buffer.data_[offset],acq->data_,sizeof(float)*2*buffer.dimensions_[0]);

					}

					printf ("%zu\n", i);

					return M;

				}


				template<class T> bool
				Write (const Matrix<T>& M, const std::string& uri) {
					return false;
				}

				template<class T> Matrix<T>
				Read (const TiXmlElement* txe) const {
					Matrix<T> M;
					return M;
				}

				template<class T> bool
				Write (const Matrix<T>& M, const TiXmlElement* txe) {
					return false;
				}

				xml_schema::properties m_props; /**< @brief Properties */


			};

		}

	};

};


#endif /* __ISMRMRD_HPP__ */
