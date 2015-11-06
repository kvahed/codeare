/*
 * HDF5File.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: kvahed
 */

#ifndef __HDF5FILE_HPP__
#define __HDF5FILE_HPP__

#include "IOFile.hpp"
#include "Matrix.hpp"
#include "Tokenizer.hpp"
#include "Workspace.hpp"

#include <boost/tokenizer.hpp>

#include <H5Cpp.h>
using namespace H5;

namespace codeare {
namespace matrix {
namespace io {

    template<class T> struct HDF5Traits;
    
    template<> struct HDF5Traits<float> {
        static PredType* PType () {
            return (PredType*) new FloatType (PredType::NATIVE_FLOAT);
        }
    };
    template<> struct HDF5Traits<double> {
        static PredType* PType () {
            return (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
        }
    };
    template<> struct HDF5Traits<cxfl> {
        static PredType* PType () {
            return (PredType*) new FloatType (PredType::NATIVE_FLOAT);
        }
    };
    template<> struct HDF5Traits<cxdb> {
        static PredType* PType () {
            return (PredType*) new FloatType (PredType::NATIVE_DOUBLE);
        }
    };
    template<> struct HDF5Traits<short> {
        static PredType* PType () {
            return (PredType*) new FloatType (PredType::NATIVE_SHORT);
        }
    };

    class HDF5File : public IOFile {

    public:

        /**
         * @brief   Open HDF5 file
         *
         * @param  fname   File name
         * @param  mode    IO mode (R/RW)
         * @param  params  Optional params
         * @param  verbose Verbosity
         */
        HDF5File  (const std::string& fname, const IOMode mode = READ,
                Params params = Params(), const bool verbose = false) :
                    IOFile(fname, mode, params, verbose), _depth(0) {
            H5::H5Library::getLibVersion(_lib_version[0],_lib_version[1],_lib_version[2]);
            Exception::dontPrint();
            try {
                m_file = H5File (fname, (mode == READ) ? H5F_ACC_RDONLY :H5F_ACC_TRUNC);
                if (this->m_verb)
                    printf ("File %s opened %s\n", fname.c_str(), (mode == READ) ? "R" : "RW");
            } catch (const FileIException& e) {
                printf ("Opening %s failed\n", fname.c_str());
                e.printError();
            }
            this->m_status = OK;

        }



        /**
         * @brief  Default destructor
         */
        virtual ~HDF5File () {
            Close ();
        }


        /**
         * @brief   Clean up and close file
         */
        virtual void Close () {
            try {
                m_file.flush(H5F_SCOPE_LOCAL);
            } catch (const Exception& e) {
                this->m_status = HDF5_ERROR_FFLUSH;
                printf ("Couldn't flush HDF5 file %s!\n%s\n", this->FileName().c_str(),
                		e.getDetailMsg().c_str());
            }
            try {
                m_file.close();
            } catch (const Exception& e) {
                this->m_status = HDF5_ERROR_FCLOSE;
                printf ("Couldn't close HDF5 file %s!\n%s\n", this->FileName().c_str(),
                		e.getDetailMsg().c_str());
            }
        }



        template<class T> Matrix<T> Read (const std::string& uri) const throw () {

            T         t       = (T) 0;
            DataSet   dataset = m_file.openDataSet(uri);
            DataSpace space   = dataset.getSpace();
            Vector<hsize_t> dims (space.getSimpleExtentNdims());
            size_t    ndim    = space.getSimpleExtentDims(&dims[0], NULL);

            if (this->m_verb) {
                printf ("Reading dataset %s ... ", uri.c_str());
                fflush(stdout);
            }

            if (is_complex(t)) {
                dims.pop_back();
                --ndim;
            }

            Vector<size_t> mdims (ndim,1);
            for (size_t i = 0; i < ndim; ++i)
                mdims[i] = dims[ndim-i-1];
            PredType* type = HDF5Traits<T>::PType();
            Matrix<T> M (mdims);
            dataset.read (&M[0], *type);

            if (this->m_verb)
                printf ("O(%s) done\n", DimsToCString(M));

            space.close();
            dataset.close();

            return M;

        }



        template<class T> bool Write (const Matrix<T>& M, const std::string& uri) throw () {

            T t = (T)0;
            Group group, *tmp;
            std::string path;

            boost::tokenizer<> tok(uri);
            std::vector<std::string> sv (Split (uri, "/"));
            std::string name = sv[sv.size() - 1];
            sv.pop_back(); // data name not part of path

            if (sv.size() == 0)
                path = "/";
            else
                for (size_t i = 0; i < sv.size(); i++) {
                    if (sv[i].compare(""))
                        path += "/";
                        path += sv[i];
                }

            if (this->m_verb)
                printf ("Creating dataset %s at path (%s)\n", name.c_str(), path.c_str());

            try {
                group = m_file.openGroup(path);
                if (this->m_verb)
                    printf ("Group %s opened for writing\n", path.c_str()) ;
            } catch (const Exception&) {
                for (size_t i = 0, depth = 0; i < sv.size(); i++) {
                    if (sv[i].compare("")) {
                        try {
                            group = (depth) ? (*tmp).openGroup(sv[i])   : m_file.openGroup(sv[i]);
                        } catch (const Exception&) {
                            group = (depth) ? (*tmp).createGroup(sv[i]) : m_file.createGroup(sv[i]);
                        }

                        tmp = &group;
                        depth++;
                    }
                }
            }

            // One more field for complex numbers
            size_t tmpdim = ndims(M);
            Vector<hsize_t> dims (tmpdim);
            for (size_t i = 0; i < tmpdim; i++)
                dims[i] = M.Dim(tmpdim-1-i);
            if (is_complex(t)) {
                dims.push_back(2);
                tmpdim++;
            }

            DataSpace space (tmpdim, &dims[0]);
            PredType*  type = HDF5Traits<T>::PType();
            DataSet set = group.createDataSet(name, (*type), space);

            set.write   (M.Ptr(), (*type));
            set.close   ();
            space.close ();

            return true;

        }


        /**
         * @brief Read a particular data set from file
         *
         * @return  Success
         */
        template<class T> Matrix<T>    Read (const TiXmlElement* txe) const {
            std::string uri (txe->Attribute("uri"));
            return this->Read<T>(uri);
        }


        /**
         * @brief  Write data to file
         *
         * @return  Success
         */
        template<class T> bool Write (const Matrix<T>& M, const TiXmlElement* txe) {
            std::string uri (txe->Attribute("uri"));
            return this->Write (M, uri);
        }

        inline void DoGroup (const H5::Group& h5g, const H5std_string& name) const {
            std::cout << std::string(_depth, ' ') << "Group: " << name << std::endl;
            _depth += 2;
            ScanAttrs(h5g);
            _depth -= 2;
        }

        inline void Read () const {
        	_depth = 4;
            std::cout << std::string(_depth, ' ') << "File name: " << m_file.getFileName() << std::endl;
            H5std_string root_str = "/";
            H5::Group h5g = m_file.openGroup(root_str);
            DoGroup(h5g,root_str);
            h5g.close();
        }

        inline void DoAttribute (const H5::Attribute& h5a) const {
            std::cout << " (" << h5a.getName();
            DataSpace space = h5a.getSpace();
            H5T_class_t dtypeclass = h5a.getTypeClass();
            hsize_t a = h5a.getStorageSize();
            switch (dtypeclass)
            {
            case H5T_INTEGER:
				break;
            case H5T_FLOAT:
            	break;
            case H5T_STRING:
            	break;
            default:
            	break;
            }
            std::cout << ")";
            space.close();
         }

        inline void ScanAttrs (const H5::Group& h5g) const {

            for (int i = 0; i < h5g.getNumAttrs(); ++i) {
                H5::Attribute h5a = h5g.openAttribute(i);
                DoAttribute(h5a);
                h5a.close();
            }

            for (hsize_t i = 0; i < h5g.getNumObjs(); ++i) {
                int otype = h5g.getObjTypeByIdx(i);
                std::string oname = h5g.getObjnameByIdx(i);
                switch(otype)
                {
                case H5G_LINK:
                    break;
                case H5G_GROUP:
                {
                    H5::Group grpid = h5g.openGroup(oname);
                    DoGroup(grpid,oname);
                    grpid.close();
                    _depth--;
                }
                break;
                case H5G_DATASET:
                {
                    H5::DataSet dsid = h5g.openDataSet(oname);
                    DoDataset(dsid,oname);
                    dsid.close();
                }
                break;
                case H5G_TYPE:
				{
					H5::DataType dtid = h5g.openDataType(oname);
					DoDatatype(dtid,oname);
					dtid.close();
				}
				break;
                default:
                    printf(" unknown?\n");
                    break;
                }
            }

        }

        inline void DoDataset (const H5::DataSet& h5d, const H5std_string& name) const {
            bool is_complex = false;
            if (_lib_version[0] >= 1 && _lib_version[1] >= 12)
            h5d.attrExists("complex");
            std::cout << std::string(_depth, ' ') << "Dataset: " << name; fflush(stdout);
            std::string wname = name;
            if (wname[0]=='/')
            	wname = wname.substr(1,wname.length());
            std::replace(wname.begin(), wname.end(), '/', '_');
            if (h5d.getTypeClass() == H5T_FLOAT) {
            	if (h5d.getFloatType() == PredType::NATIVE_FLOAT) {
            		if (is_complex) {
            			Matrix<cxfl> M = Read<cxfl>(name);
            			wspace.Add(wname,M);
            		} else {
            			Matrix<float> M = Read<float>(name);
            			wspace.Add(wname,M);
            		}
            	} else if (h5d.getFloatType() == PredType::NATIVE_DOUBLE) {
            		if (is_complex) {
            			Matrix<cxdb> M = Read<cxdb>(name);
            			wspace.Add(wname,M);
            		} else {
            			Matrix<double> M = Read<double>(name);
            			wspace.Add(wname,M);
            		}
            	}
            }
            std::cout << std::endl;
        }


        inline void DoDatatype (const H5::DataType& h5t, const H5std_string& name) const {}

    private:

        HDF5File (const HDF5File&) : _depth(0) {}
        HDF5File  () : _depth(0) {}
        H5File m_file; /// @brief My file
        mutable size_t _depth;
        unsigned lib_version[3];


    };

}// namespace io
}// namespace matrix
}// namespace codeare



    template<class T> inline static bool _h5write (const Matrix<T>& M, const std::string& fname,
    		const std::string& uri) {
        using namespace codeare::matrix::io;
        HDF5File h5f (fname, WRITE);
        h5f.Write(M, uri);
        return true;
    }
#define h5write(X,Y) _h5write (X,Y,#X)

    template<class T> inline static Matrix<T> h5read (const std::string& fname,
    		const std::string& uri) {
        using namespace codeare::matrix::io;
        HDF5File h5f (fname, READ);
        return h5f.Read<T>(uri);
    }

#endif /* __HDF5FILE_HPP__ */
