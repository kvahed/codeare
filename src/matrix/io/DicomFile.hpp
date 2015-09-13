/*
 * Dicom.hpp
 *
 *  Created on: Sep 11, 2015
 *      Author: kvahed
 */

#ifndef SRC_MATRIX_IO_DICOM_HPP_
#define SRC_MATRIX_IO_DICOM_HPP_

#include "IOFile.hpp"

#include "itkVersion.h"

#include "itkImage.h"
#include "itkMinimumMaximumImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"

#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"

#include <itksys/SystemTools.hxx>

#if ITK_VERSION_MAJOR >= 4
#include "gdcmUIDGenerator.h"
#else
#include "gdcm/src/gdcmFile.h"
#include "gdcm/src/gdcmUtil.h"
#endif

namespace codeare {
namespace matrix {
namespace io {

class DicomFile : public IOFile {
public:

    typedef uint16_t PixelType;
	typedef itk::NumericSeriesFileNames OutputNamesGeneratorType;
    typedef itk::GDCMImageIO ImageIOType;

	DicomFile  (const std::string& fname, const IOMode mode = READ,
			Params params = Params(), const bool verbose = false) :
				IOFile(fname, mode, params, verbose) {

        // Directory create / read / write success
		std::stringstream testfilestr;
		testfilestr << std::string(".");
		testfilestr << fname;
		ofstream testfile(testfilestr.str());
		if (!testfile.is_open())
			throw OPEN_RW_FAILED;
		std::remove(testfilestr.str().c_str());

		// Create outputfile generator
        if (mode == WRITE) {
            _ofnames = OutputNamesGeneratorType::New();
        }
        _gdcmIO = ImageIOType::New();
        
		this->m_status = OK;
        
	}

	template<class T> Matrix<T> Read (const std::string& uri) throw () {
		return Matrix<T>();
	}

	template<class T> bool Write (const Matrix<T>& M, const std::string& uri) throw () {
        
        typedef itk::Image<PixelType, 2> InputImageType;
        typedef itk::Image<PixelType, 2> OutputImageType;
        typedef itk::ImageSeriesWriter< InputImageType, OutputImageType > SeriesWriterType;
        
        // Make the output directory and generate the file names.
        itksys::SystemTools::MakeDirectory(uri);
        size_t n_images = numel(M)/(size(M,0)*size(M,1));

        // Generate the file names
        std::string seriesFormat("");
        seriesFormat = seriesFormat + "/" + "IM%d.dcm";
        _ofnames->SetSeriesFormat (seriesFormat.c_str());
        _ofnames->SetStartIndex (1);
        _ofnames->SetEndIndex (numel(M)/(size(M,0)*size(M,1)));

        // Prepare writer for this series

        SeriesWriterType::Pointer _writer = SeriesWriterType::New();
        _writer->SetImageIO( _gdcmIO );
        _writer->SetFileNames( _ofnames->GetFileNames() );

        // Write
        try {
            _writer->Update();
        } catch (const itk::ExceptionObject & excp ) {
            std::cerr << "Exception thrown while writing the series " << std::endl;
            std::cerr << excp.what() << std::endl;
            return false;
        }
        
        std::cout << "  Wrote " << n_images << "dicom images to " << uri.c_str();
        return true;
        
    }
    
    
	/**
	 * @brief  Default destructor
	 */
	virtual ~DicomFile () {}

private:
	OutputNamesGeneratorType::Pointer _ofnames;
    ImageIOType::Pointer _gdcmIO;
};

} /* namespace io */
} /* namespace matrix */
} /* namespace codeare */

#endif /* SRC_MATRIX_IO_DICOM_HPP_ */
