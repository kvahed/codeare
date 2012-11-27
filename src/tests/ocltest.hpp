/**************
 ** includes **
 **************/
 
// CoDEARE
# include "matrix/ocl/oclmatrix/oclMatrix.hpp"
# include "matrix/ocl/oclConnection.hpp"
# include "matrix/ocl/oclmatrix/oclIO.hpp"

// own
# include "matrix/ocl/timer.h"


/**********************
 ** type definitions **
 **********************/
typedef float elem_type;
typedef oclMatrix <elem_type> Matrix_type;




template <class T>
bool
mat_equal             ( const oclMatrix <T> & ocl_mat,
                        const    Matrix <T> &    smat )
{

  ocl_mat.getData (); // !!! //
  Matrix <bool> mat_comp = (smat == ocl_mat);
  bool result = true;
  for (int i = 0; i < smat.Height (); i++)
    for (int j = 0; j < smat.Width (); j++)
      result &= mat_comp (i,j);
  return result;

}


template <class T>
bool
oclCreatorsTest       ( bool verbose )
{

  std::cout <<                           std::endl;
  std::cout << " ****************** " << std::endl;
  std::cout << " ** constructors ** " << std::endl;
  std::cout << " ****************** " << std::endl;
  std::cout <<                           std::endl;

  /* default constructor */
  {
    if (verbose) std::cout << " * oclMatrix ( )                       ";
    oclMatrix <T> ocl_mat;
       Matrix <T>     mat;
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with dimensions */
  {
    if (verbose) std::cout << " * oclMatrix ( dims )                  ";
    size_t dims [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
       Matrix <T>     mat (dims);
    oclMatrix <T> ocl_mat (dims);
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with dimensions and resolutions */
  {
    if (verbose) std::cout << " * oclMatrix ( dims, res )             ";
    size_t dims [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
     float  res [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
       Matrix <T>     mat (dims, res);
    oclMatrix <T> ocl_mat (dims, res);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct sqare with dimension */
  {
    if (verbose) std::cout << " * oclMatrix ( n )                     ";
    size_t n = 16;
       Matrix <T>     mat (n);
    oclMatrix <T> ocl_mat (n);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }
  
  /* construct with rows and columns */
  {
    if (verbose) std::cout << " * oclMatrix ( rows, cols )            ";
    size_t rows = 16,
           cols = 16;
       Matrix <T>     mat (rows, cols);
    oclMatrix <T> ocl_mat (rows, cols);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with rows, columns and slices */
  {
    if (verbose) std::cout << " * oclMatrix ( rows, cols, slices )    ";
    size_t rows   = 16,
           cols   = 16,
           slices = 16;
       Matrix <T>     mat (rows, cols, slices);
    oclMatrix <T> ocl_mat (rows, cols, slices);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with all dimensions */
  {
    if (verbose) std::cout << " * oclMatrix ( col, lin, cha, ... )    ";
       Matrix <T>     mat (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    oclMatrix <T> ocl_mat (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }
  
  /* construct with zeros */
  {
    if (verbose) std::cout << " * zeros ( col, lin, cha, ... )        ";
       Matrix <T>     mat = zeros <T> (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    oclMatrix <T> ocl_mat = zeros <T> (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    assert (mat_equal <T> (ocl_mat, mat) && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

}


template <class T>
oclMatrix <T>
oclInit2D              ( size_t dimX, size_t dimY, const T & offset = 0 )
{

  oclMatrix <T> ocl_mat (dimX, dimY);
  for (size_t i = 0; i < dimX; i++)
    for (size_t j = 0; j < dimY; j++)
      ocl_mat (i, j) = i + j + offset;
      
  return ocl_mat;

}


template <class T>
Matrix <T>
Init2D              ( size_t dimX, size_t dimY, const T & offset = 0 )
{

  Matrix <T> mat (dimX, dimY);
  for (size_t i = 0; i < dimX; i++)
    for (size_t j = 0; j < dimY; j++)
      mat (i, j) = i + j + offset;
      
  return mat;

}


template <class T>
bool
oclArithmeticsTest    ( bool verbose )
{

  std::cout <<                          std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout << " ** arithmetics ** " << std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout <<                          std::endl;

  size_t  dimX     = 1024,
          dimY     = 1024;
  double  time_ocl = 0.0,
          time_s   = 0.0;
  MyTimer t;  
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl;
  
  
  /* add two matrices */
  {
    if (verbose) std::cout << " * res = m1 + m2               ";
    t.tic (time_default);
      Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
      Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
      Matrix <T>     res                  (dimX, dimY);
      res = mat1 + mat2;
    time_s = t.tic (time_default);
      oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
      oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
      oclMatrix <T> ocl_res                  (dimX, dimY);
      ocl_res = ocl_mat1 + ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* subtract two matrices */
  {
    if (verbose) std::cout << " * res = m1 - m2               ";
    t.tic (time_default);
      Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
      Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
      Matrix <T>     res                  (dimX, dimY);
      res = mat1 - mat2;
    time_s = t.tic (time_default);
      oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
      oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
      oclMatrix <T> ocl_res                  (dimX, dimY);
      ocl_res = ocl_mat1 - ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* increment matrix (uniform) */
  {
    if (verbose) std::cout << " * m += scalar                 ";
    const T scalar = (T) 33.3;
    t.tic (time_default);
      Matrix <T>     mat    =    Init2D <T> (dimX, dimY);
      mat     += scalar;
    time_s = t.tic (time_default);
      oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
      ocl_mat += scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* decrement matrix (uniform) */
  {
    if (verbose) std::cout << " * m -= scalar                 ";
    const T scalar = (T) 33.3;
    t.tic (time_default);
      Matrix <T>     mat    =    Init2D <T> (dimX, dimY);
      mat     -= scalar;
    time_s   = t.tic (time_default);
      oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
      ocl_mat -= scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* increment matrix (non uniform) */
  {
    if (verbose) std::cout << " * m += matrix                 ";
    const    Matrix <T>     mat_inc =    Init2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_inc = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <T>        mat =    Init2D <T> (dimX, dimY);
      mat     +=     mat_inc;
    time_s   = t.tic (time_default);
      oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
      ocl_mat += ocl_mat_inc;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* decrement matrix (non uniform) */
  {
    if (verbose) std::cout << " * m -= matrix                 ";
    const    Matrix <T>     mat_dec =    Init2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_dec = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <T>        mat =    Init2D <T> (dimX, dimY);
      mat     -=     mat_dec;
    time_s   = t.tic (time_default);
      oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
      ocl_mat -= ocl_mat_dec;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

}




/*******************
 ** test function **
 *******************/
template <class T> bool
oclmatrixtest (Connector<T>* rc) {

  oclCreatorsTest <elem_type>    (true);
  oclArithmeticsTest <elem_type> (true);

  // choose tests
  enum Test_Type {constructors, ocl_scalar_add, ocl_mat_inc, ocl_mat_add, ocl_mat_sub, ocl_scalar_sub};
  bool tests [6] = {      true,           true,        true,        true,        true,           true};
  
  
  std::cout << std::endl;
  std::cout << " * oclmatrixtest " << std::endl << std::endl;
  
  // choose dimensions
  int dimX = 512, dimY = 512;

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
  
  if (tests [ocl_scalar_add])
  {
    std::cout << " * scalar addition" << std::endl;
    mat_zeros += 6;
  }
  
  if (tests [ocl_scalar_sub])
  {
    std::cout << " * scalar subraction" << std::endl;
    mat_zeros -= 6;
  }

  if (tests [ocl_mat_inc])
  {
    std::cout << " * oclMatrix <T> addition" << std::endl;
    mat_zeros += mat_zeros;
  }
  
///////////////////////////////////////
  
  MyTimer t;
  double time_ocl, time_s;
  Matrix <elem_type> smat;
  Matrix <bool> mat_comp;
  bool result;
  
  mat_zeros.getData (); 
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
 
  /* print list of loaded objects */
  oclConnection :: Instance () -> print_ocl_objects ();
  
  std::cout << " * ocltest finished!" << std::endl;

	return true;

}
