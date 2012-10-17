/**************
 ** includes **
 **************/
 
// CoDEARE
# include "matrix/ocl/oclmatrix/oclMatrix.hpp"
# include "matrix/ocl/oclmatrix/oclCreators.hpp"
# include "matrix/Creators.hpp"



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


  std::cout << std::endl;
  std::cout << " * oclmatrixtest " << std::endl << std::endl;
  
  int dimX = 10, dimY = 1;

  // constructors
  Matrix_type mat;
  Matrix <size_t> sz (2, 1);
  sz [0] = dimX, sz[1] = dimY;
  Matrix_type mat_zeros = zeros <elem_type> (Matrix <size_t> (sz));
  
  Matrix_type mat_copy (mat_zeros);
  std::cout << " print \"mat_copy\": " << std::endl << mat_copy << std::endl;
  
  mat_zeros += 3;
  
  mat_zeros += (Matrix <elem_type>) mat_zeros;
  
  mat_zeros = mat_zeros + mat_zeros;
  
  std::cout << " print \"mat_zeros\": " << std::endl << mat_zeros << std::endl;
  
  std::cout << " * finished!" << std::endl;

	return true;

}
