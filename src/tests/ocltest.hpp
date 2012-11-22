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
  bool tests [5] = {      true,         true,     false,        true,        true};
  
  
  std::cout << std::endl;
  std::cout << " * oclmatrixtest " << std::endl << std::endl;
  
  // choose dimensions
  int dimX = 10000, dimY = 2000;

  bool verbose;
  if (dimX + dimY < 100)
    verbose = true;

  Matrix_type mat_zeros;
  if (verbose)
    std::cout << "mat_zeros" << mat_zeros << std::endl;
  mat_zeros.getData ();

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
//    std::cout << " print \"mat_copy\": " << std::endl << mat_copy << std::endl;
  }
  
  if (tests [m_scalar_add])
  {
    std::cout << " * scalar addition" << std::endl;
    mat_zeros.getData (); /* !!! */
    mat_zeros += 3;
  }
/*
  Matrix_type mat;
  std::cout << " !!!!!!! mat_zeros:" << std::endl;
  std::cout << mat_zeros << std::endl;
  mat_zeros + mat_zeros;
  mat = mat_zeros;
  std::cout << " !!!!!!! mat:" << std::endl;
  std::cout << mat << std::endl;
*/  
/*  if (tests [m_mat_add])
  {
    std::cout << " * Matrix<T> addition" << std::endl;
    mat_zeros += (Matrix <elem_type>) mat_zeros;
  }
*/  
///////////////////////////////////////
  
  MyTimer t;
  double time_ocl, time_s;
  Matrix <elem_type> smat;
  Matrix <bool> mat_comp;
  bool result;
  
  smat = mat_zeros;
  
  if (tests [ocl_mat_add])
  {
    
    Matrix_type mat1 = mat_zeros, mat2 = mat_zeros;
    if (verbose)
      std::cout << "mat_zeros\n" << mat_zeros << std::endl;
    //////
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << " * oclMatrix<T> addition" << std::endl;

    t.tic (time_default);
    mat_zeros = mat1 + mat2;
    time_ocl = t.tic (time_default);

    if (verbose)
      std::cout << " print \"mat_zeros\": " << std::endl << mat_zeros << std::endl;

    t.tic (time_default);
    smat = smat + smat;
    time_s = t.tic (time_default);

    mat_zeros.getData (); // !!! //
    Matrix <bool> mat_comp = (smat == mat_zeros);
    result = true;
    for (int i = 0; i < dimX; i++)
      for (int j = 0; j < dimY; j++)
        result &= mat_comp (i,j);
    assert (result == true);
    std::cout << "time_ocl: " << time_ocl << ", time_s: " << time_s << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
  }

  ///////

  if (tests [ocl_mat_sub])
  {
  
    std::cout << " * oclMatrix<T> subtraction" << std::endl;

    t.tic (time_default);
    mat_zeros = mat_zeros - mat_zeros;
    time_ocl = t.tic (time_default);
  
    if (verbose)
      std::cout << " print \"mat_zeros\": " << std::endl << mat_zeros << std::endl;
  
    t.tic (time_default);
    smat = smat - smat;
    time_s = t.tic (time_default);
  
    mat_zeros.getData (); // !!! //
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
 /*
 
  // choose dimensions
  int dimX = 2, dimY = 2;
  
  Matrix <size_t> sz (2, 1);
  sz [0] = dimX, sz[1] = dimY;
  Matrix_type mat1 = zeros <elem_type> (sz);
  Matrix_type mat2 = zeros <elem_type> (sz);
  
  mat1 += 3;
  mat2 += 2;
  
  std::cout << "mat1: \n" << mat1 << std::endl;
  std::cout << "mat2: \n" << mat2 << std::endl;
  
  Matrix_type result = mat1 + mat2;
  
  result += 0.5;

  std::cout << "mat1: \n" << mat1 << std::endl;
  std::cout << "mat2: \n" << mat2 << std::endl;
  std::cout << "result: \n" << result << std::endl;
    
  mat1 = result + mat2;
  
  std::cout << "mat1: \n" << mat1 << std::endl;
  std::cout << "mat2: \n" << mat2 << std::endl;
  std::cout << "result: \n" << result << std::endl;
*/
 
  oclConnection :: Instance () -> print_ocl_objects ();
  
  std::cout << " * finished!" << std::endl;

	return true;

}
