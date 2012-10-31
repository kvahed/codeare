/**************
 ** includes **
 **************/
 
// CoDEARE
# include "matrix/ocl/oclmatrix/oclMatrix.hpp"
# include "matrix/ocl/oclConnection.hpp"
# include "matrix/ocl/oclmatrix/oclIO.hpp"
# include "matrix/Creators.hpp"

// own
# include "matrix/ocl/timer.h"


/**********************
 ** type definitions **
 **********************/
typedef float elem_type;
typedef oclMatrix <elem_type> Matrix_type;



/*******************
 ** test function **
 *******************/
template <class T> bool
oclmatrixtest (Connector<T>* rc) {

  // choose tests
  enum Test_Type {constructors, m_scalar_add, m_mat_add, ocl_mat_add, ocl_mat_sub};
  bool tests [5] = {      true,        true,      true,        true,       false};
  
  
  std::cout << std::endl;
  std::cout << " * oclmatrixtest " << std::endl << std::endl;
  
  // choose dimensions
  int dimX = 10, dimY = 10;

  Matrix_type mat_zeros;
  std::cout << "mat_zeros" << std::endl;


  if (tests [constructors])
  {
    // constructors
    Matrix_type mat;
    std::cout << "mat" << std::endl;;
    Matrix <size_t> sz (2, 1);
    sz [0] = dimX, sz[1] = dimY;
    mat_zeros = zeros <elem_type> (Matrix <size_t> (sz));
  
    Matrix_type mat_copy (mat_zeros);
    std::cout << "mat_copy" << std::endl;
    std::cout << " print \"mat_copy\": " << std::endl << mat_copy << std::endl;
  }
  
  if (tests [m_scalar_add])
  {
    std::cout << " * scalar addition" << std::endl;
    mat_zeros += 3;
  }
  
  if (tests [m_mat_add])
  {
    std::cout << " * Matrix<T> addition" << std::endl;
    mat_zeros += (Matrix <elem_type>) mat_zeros;
  }
  
///////////////////////////////////////
  
  MyTimer t;
  double time_ocl, time_s;
  Matrix <elem_type> smat;
  Matrix <bool> mat_comp;
  bool result;
  
  if (tests [ocl_mat_add])
  {
    smat = mat_zeros;
    //////
    std::cout << " * oclMatrix<T> addition" << std::endl;

    t.tic (time_default);
    mat_zeros = mat_zeros + mat_zeros;
    time_ocl = t.tic (time_default);

    std::cout << " print \"mat_zeros\": " << std::endl << mat_zeros << std::endl;

    t.tic (time_default);
    smat = smat + smat;
    time_s = t.tic (time_default);
    
    Matrix <bool> mat_comp = (smat == mat_zeros);
    result = true;
    for (int i = 0; i < dimX; i++)
      for (int j = 0; j < dimY; j++)
        result &= mat_comp (i,j);
    assert (result == true);
    std::cout << "time_ocl: " << time_ocl << ", time_s: " << time_s << std::endl;
  }

  ///////

  if (tests [ocl_mat_sub])
  {
  
    std::cout << " * oclMatrix<T> subtraction" << std::endl;

    t.tic (time_default);
    mat_zeros = mat_zeros - mat_zeros;
    time_ocl = t.tic (time_default);
  
    t.tic (time_default);
    smat = smat - smat;
    time_s = t.tic (time_default);
  
    mat_comp = (smat == mat_zeros);
    result = true;
    for (int i = 0; i < dimX; i++)
      for (int j = 0; j < dimY; j++)
        result &= mat_comp (i,j);
    assert (result == true);
    std::cout << "time_ocl: " << time_ocl << ", time_s: " << time_s << std::endl;
  }
  
//////////////////////////////////////////
  
  
//  std::cout << " print \"mat_zeros\": " << std::endl << mat_zeros << std::endl;
  
  std::cout << " * finished!" << std::endl;

	return true;

}
