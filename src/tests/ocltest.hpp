/**************
 ** includes **
 **************/
 
// CoDEARE
# include "matrix/ocl/oclmatrix/oclMatrix.hpp"
# include "matrix/ocl/oclConnection.hpp"
# include "matrix/ocl/oclmatrix/oclIO.hpp"

// C++ std lib
# include <iomanip>

// own
# include "matrix/ocl/timer.h"


/**********************
 ** type definitions **
 **********************/
typedef double elem_type;
typedef oclMatrix <elem_type> Matrix_type;


/*************
 ** structs **
 *************/

struct result
{
  
  bool equal;
  elem_type mean_abs_err;
  
};


template <class T>
result
mat_equal             ( const oclMatrix <T> & ocl_mat,
                        const    Matrix <T> &    smat )
{

  ocl_mat.getData (); // !!! //
  Matrix <bool> mat_comp = (smat == ocl_mat);
  result res = {true, 0.0};
  for (int i = 0; i < smat.Height (); i++)
    for (int j = 0; j < smat.Width (); j++)
    {
      res.equal &= mat_comp (i,j);
    }
  
  if (! res.equal)
  {
    Matrix <elem_type> diff = smat - ((Matrix <elem_type>) ocl_mat);
    for (int i = 0; i < smat.Height (); i++)
      for (int j = 0; j < smat.Width (); j++)
        res.mean_abs_err += abs (diff (i, j));
    res.mean_abs_err /= smat.Size ();
  }
  
  return res;

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
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with dimensions */
  {
    if (verbose) std::cout << " * oclMatrix ( dims )                  ";
    size_t dims [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
       Matrix <T>     mat (dims);
    oclMatrix <T> ocl_mat (dims);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with dimensions and resolutions */
  {
    if (verbose) std::cout << " * oclMatrix ( dims, res )             ";
    size_t dims [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
     float  res [16] = {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2};
       Matrix <T>     mat (dims, res);
    oclMatrix <T> ocl_mat (dims, res);
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct sqare with dimension */
  {
    if (verbose) std::cout << " * oclMatrix ( n )                     ";
    size_t n = 16;
       Matrix <T>     mat (n);
    oclMatrix <T> ocl_mat (n);
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }
  
  /* construct with rows and columns */
  {
    if (verbose) std::cout << " * oclMatrix ( rows, cols )            ";
    size_t rows = 16,
           cols = 16;
       Matrix <T>     mat (rows, cols);
    oclMatrix <T> ocl_mat (rows, cols);
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
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
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

  /* construct with all dimensions */
  {
    if (verbose) std::cout << " * oclMatrix ( col, lin, cha, ... )    ";
       Matrix <T>     mat (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    oclMatrix <T> ocl_mat (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }
  
  /* construct with zeros */
  {
    if (verbose) std::cout << " * zeros ( col, lin, cha, ... )        ";
       Matrix <T>     mat = zeros <T> (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    oclMatrix <T> ocl_mat = zeros <T> (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2);
    assert (mat_equal <T> (ocl_mat, mat).equal && "Test failed! ");
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
      ocl_mat (i, j) = (i + j + offset) / ((dimX + dimY) / 2);
      
  return ocl_mat;

}


template <class T>
Matrix <T>
Init2D              ( size_t dimX, size_t dimY, const T & offset = 0 )
{

  Matrix <T> mat (dimX, dimY);
  for (size_t i = 0; i < dimX; i++)
    for (size_t j = 0; j < dimY; j++)
      mat (i, j) = (i + j + offset) / ((dimX + dimY) / 2);
      
  return mat;

}


template <class T>
bool
oclAccessTest       ( bool verbose )
{

  std::cout <<                     std::endl;
  std::cout << " ************ " << std::endl;
  std::cout << " ** access ** " << std::endl;
  std::cout << " ************ " << std::endl;
  std::cout <<                     std::endl;
  
  size_t dimX = 512,
         dimY = 512;

  /* elementwise access (rvalue) */
  {
    if (verbose) std::cout << " * At (p), [p] (rvalue)       ";
       Matrix <T>    mat =    Init2D <T> (dimX, dimY);
    oclMatrix <T> oclmat = oclInit2D <T> (dimX, dimY);
    bool test = true;
    for (int i = 0; i < mat.Size (); i++)
      test &= (mat.At (i) == oclmat.At (i)) && (mat [i] == oclmat [i]);
    assert (test && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }
  

  /* elementwise access (lvalue) */
  {
    if (verbose) std::cout << " * At (p), [p] (lvalue)       ";
       Matrix <T>    mat =    Init2D <T> (dimX, dimY);
    oclMatrix <T> oclmat = oclInit2D <T> (dimX, dimY);
    bool test = true;
    for (int i = 0; i < mat.Size (); i++)
    {
         mat.At (i) =    mat [i] = i;
      oclmat.At (i) = oclmat [i] = i;
      test &= (mat.At (i) == oclmat.At (i)) && (mat [i] == oclmat [i]);
    }
    assert (test && " Test failed! ");
    if (verbose) std::cout << "passed!" << std::endl;
  }

}


template <class T>
bool
oclComparisonsTest    ( bool verbose )
{

  std::cout <<                          std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout << " ** comparisons ** " << std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout <<                          std::endl;

  size_t  dimX     = 2048,
          dimY     = 2048;
  double  time_ocl = 0.0,
          time_s   = 0.0;
  MyTimer t;  
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl;

  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A == s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <bool>   res = mat1 == scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 == scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A != s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <bool>   res = mat1 != scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 != scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A >  s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <bool>   res = mat1 > scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 > scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }
  
  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A >= s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <bool>   res = mat1 >= scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 >= scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A <  s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <bool>   res = mat1 < scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 < scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }
  
  /* compare matrix and scalar */
  {
    if (verbose) std::cout << " * A <= s                      ";
           T       scalar = 3;
    const    Matrix <T>        mat1 = Init2D <T> (dimX, dimY);
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
         Matrix <bool>      res =     mat1 <= scalar;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 <= scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrices */
  {
    if (verbose) std::cout << " * A == B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                            mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 == mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 == ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrices */
  {
    if (verbose) std::cout << " * A != B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 != mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 != ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }
  
  /* compare matrices */
  {
    if (verbose) std::cout << " * A >  B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 > mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 > ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrices */
  {
    if (verbose) std::cout << " * A >= B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 >= mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 >= ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }
  
  /* compare matrices */
  {
    if (verbose) std::cout << " * A <  B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 < mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 < ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* compare matrices */
  {
    if (verbose) std::cout << " * A <= B                      ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 <= mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 <= ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise AND */
  {
    if (verbose) std::cout << " * C = A && B                  ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 && mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 && ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise OR */
  {
    if (verbose) std::cout << " * C = A || B                  ";
    int pos_x = 10,
        pos_y = 11;
    const    Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
             Matrix <T>     mat2 = Init2D <T> (dimX, dimY);
                          mat2 (pos_x, pos_y) += 1;
    const oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T>    ocl_mat2 = oclInit2D <T> (dimX, dimY);
                           ocl_mat2 (pos_x, pos_y) += 1;
    t.tic (time_default);
      Matrix <bool>   res = mat1 || mat2;
    time_s = t.tic (time_default);
      oclMatrix <bool>  ocl_res = ocl_mat1 || ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <bool> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

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

  size_t  dimX     = 2048,
          dimY     = 2048;
  double  time_ocl = 0.0,
          time_s   = 0.0;
  MyTimer t;  
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl;
  
  
 // verbosity = VERB_HIGH;


  /* add two matrices */
  {
    if (verbose) std::cout << " * C = A + B                   ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
          res = mat1 + mat2;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 + ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;   
  }

  /* subtract two matrices */
  {
    if (verbose) std::cout << " * C = A - B                   ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = mat1 - mat2;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 - ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* increment matrix (uniform) */
  {
    if (verbose) std::cout << " * M = A + scalar              ";
    const T scalar = (T) 33.3;
       Matrix <T>     mat,
                      mat_A =    Init2D <T> (dimX, dimY);
    oclMatrix <T> ocl_mat,
                  ocl_mat_A = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     =     mat_A + scalar;
    time_s = t.tic (time_default);
      ocl_mat = ocl_mat_A + scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* decrement matrix (uniform) */
  {
    if (verbose) std::cout << " * M = A - scalar              ";
    const T scalar = (T) 33.3;
       Matrix <T>     mat,
                      mat_A =    Init2D <T> (dimX, dimY);
    oclMatrix <T> ocl_mat,
                  ocl_mat_A = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     =     mat_A + scalar;
    time_s = t.tic (time_default);
      ocl_mat = ocl_mat_A + scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* increment matrix (uniform) */
  {
    if (verbose) std::cout << " * M += scalar                 ";
    const T scalar = (T) 33.3;
       Matrix <T>     mat =    Init2D <T> (dimX, dimY);
    oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     += scalar;
    time_s = t.tic (time_default);
      ocl_mat += scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* decrement matrix (uniform) */
  {
    if (verbose) std::cout << " * M -= scalar                 ";
    const T scalar = (T) 33.3;
       Matrix <T>     mat =    Init2D <T> (dimX, dimY);
    oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     -= scalar;
    time_s   = t.tic (time_default);
      ocl_mat -= scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* increment matrix (non uniform) */
  {
    if (verbose) std::cout << " * A += B                      ";
             Matrix <T>         mat =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat_inc =    Init2D <T> (dimX, dimY);
          oclMatrix <T>     ocl_mat = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_inc = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     +=     mat_inc;
    time_s   = t.tic (time_default);
      ocl_mat += ocl_mat_inc;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* decrement matrix (non uniform) */
  {
    if (verbose) std::cout << " * A -= B                      ";
             Matrix <T>         mat =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat_dec =    Init2D <T> (dimX, dimY);
          oclMatrix <T>     ocl_mat = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_dec = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      mat     -=     mat_dec;
    time_s   = t.tic (time_default);
      ocl_mat -= ocl_mat_dec;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_mat, mat).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* raise to higher power */
  {
    if (verbose) std::cout << " * C = A ^ p                   ";
    const T p = 1.55;
    const    Matrix <T>     mat_A =    Init2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_A = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      Matrix <T>        mat_C =     mat_A ^ p;
    time_s   = t.tic (time_default);
      oclMatrix <T> ocl_mat_C = ocl_mat_A ^ p;
    time_ocl = t.tic (time_default);
    result res = mat_equal <T> (ocl_mat_C, mat_C);
    assert ((res.equal || res.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << res.mean_abs_err << fixed << std::endl;
  }

  /* raise to higher power */
  {
    if (verbose) std::cout << " * M ^= p                      ";
    const T p = 1.55;
       Matrix <T>     mat =    Init2D <T> (dimX, dimY);
    oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
          mat ^= p;
    time_s   = t.tic (time_default);
      ocl_mat ^= p;
    time_ocl = t.tic (time_default);
    result res = mat_equal <T> (ocl_mat, mat);
    assert ((res.equal || res.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << res.mean_abs_err << fixed << std::endl;
  }

  /* elementwise multiplication of matrix and scalar */
  {
    if (verbose) std::cout << " * C = A * s                   ";
    const     T       scalar = 3.33;
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = mat1 * scalar;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 * scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise multiplacation of matrix and scalar */
  {
    if (verbose) std::cout << " * C *= s                      ";
    const     float       scalar = 3;
             Matrix <T>     res   =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res   = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      res *= scalar;
    time_s = t.tic (time_default);
      ocl_res *= scalar;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise multiplication of two matrices */
  {
    if (verbose) std::cout << " * C = A .* B                  ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = mat1 * mat2;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 * ocl_mat2;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise multiplacation of two matrices */
  {
    if (verbose) std::cout << " * C .*= A                     ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res  =    Init2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res  = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      res *= mat1;
    time_s = t.tic (time_default);
      ocl_res *= ocl_mat1;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }

  /* elementwise division of matrix by scalar */
  {
    if (verbose) std::cout << " * C = A / s                   ";
    const     T       scalar = 3.33;
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = mat1 / scalar;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 / scalar;
    time_ocl = t.tic (time_default);
    result r = mat_equal <T> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed << std::endl;
  }

  /* elementwise division of matrix by scalar */
  {
    if (verbose) std::cout << " * C /= s                      ";
    const     T      scalar = 3.33;
             Matrix <T>     res  =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res  = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      res /= scalar;
    time_s = t.tic (time_default);
      ocl_res /= scalar;
    time_ocl = t.tic (time_default);
    result r = mat_equal <T> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed << std::endl;
  }

  /* elementwise division of two matrices */
  {
    if (verbose) std::cout << " * C = A ./ B                  ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat2 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat2 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = mat1 / mat2;
    time_s = t.tic (time_default);
      ocl_res = ocl_mat1 / ocl_mat2;
    time_ocl = t.tic (time_default);
    result r = mat_equal <T> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 1e-0) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed << std::endl;
  }

  /* elementwise division of two matrices */
  {
    if (verbose) std::cout << " * C ./= A                     ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res  =    Init2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res  = oclInit2D <T> (dimX, dimY);
    t.tic (time_default);
      res = res / mat1;
    time_s = t.tic (time_default);
      ocl_res /= ocl_mat1;
    time_ocl = t.tic (time_default);
    result r = mat_equal <T> (ocl_res, res);
    assert ((r.equal || r.mean_abs_err < 1e-0) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << r.mean_abs_err << fixed << std::endl;
  }
  
  /* unary minus */
  {
    if (verbose) std::cout << " * C = -A                      ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = -mat1;
    time_s = t.tic (time_default);
      ocl_res = -ocl_mat1;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

  /* unary plus */
/*  {
    if (verbose) std::cout << " * C = +A                     ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
      res = +mat1;
    time_s = t.tic (time_default);
      ocl_res = +ocl_mat1;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;    
  }*/

}



template <class T>
bool
oclLinalgTest           ( bool verbose )
{

  std::cout <<                          std::endl;
  std::cout << " ******************** " << std::endl;
  std::cout << " ** Linear Algebra ** " << std::endl;
  std::cout << " ******************** " << std::endl;
  std::cout <<                          std::endl;

  size_t  dimX     = 2048,
          dimY     = 2048;
  double  time_ocl = 0.0,
          time_s   = 0.0;
  MyTimer t;  
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl;

  /* transpose matrix */
  {
    if (verbose) std::cout << " * C = ! A                     ";
    const    Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
             Matrix <T>     res                  (dimX, dimY);
    const oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
          oclMatrix <T> ocl_res                  (dimX, dimY);
    t.tic (time_default);
          res = ! mat1;
    time_s = t.tic (time_default);
      ocl_res = ! ocl_mat1;
    time_ocl = t.tic (time_default);
    assert (mat_equal <T> (ocl_res, res).equal && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;   
  }

  /* matrix product */
  {
    if (verbose) std::cout << " * C = A * B                   ";
    const    Matrix <T>     mat_A =    Init2D <T> (dimX, dimY);
    const    Matrix <T>     mat_B =    Init2D <T> (dimY, dimX);
    const oclMatrix <T> ocl_mat_A = oclInit2D <T> (dimX, dimY);
    const oclMatrix <T> ocl_mat_B = oclInit2D <T> (dimY, dimX);        
    t.tic (time_default);
      Matrix <T>        mat_C =     mat_A ->*     mat_B;
    time_s   = t.tic (time_default);
      oclMatrix <T> ocl_mat_C = ocl_mat_A ->* ocl_mat_B;
    time_ocl = t.tic (time_default);
    result res = mat_equal <T> (ocl_mat_C, mat_C);
    assert ((res.equal || res.mean_abs_err < 1e-1) && " Test failed! ");
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)"
                           << " -> mean abs. err. = " << scientific << res.mean_abs_err << fixed << std::endl;
  }

}



template <class T>
void
exec_tests      (bool verbose)
{

  oclCreatorsTest     <T> (verbose);
  oclArithmeticsTest  <T> (verbose);
  oclLinalgTest       <T> (verbose);
  oclAccessTest       <T> (verbose);
  oclComparisonsTest  <T> (verbose);

}



/*******************
 ** test function **
 *******************/
template <class T> bool
oclmatrixtest (Connector<T>* rc) {

  /* format output */
  std::cout << fixed;
  std::cout << setprecision (2);

  try
  {
    std::cout << std::endl << std::endl;
    std::cout << " ----------- // float  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <float>  (true);
    std::cout << std::endl << std::endl;
    std::cout << " ----------- // double \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <double> (true);
    std::cout << std::endl;
  }
  catch (oclError & err)
  {
    print_optional (err, VERB_NONE);
  }
 
  /* print list of loaded objects */
  std::cout << std::endl;
  std::cout << " * Data objects in memory: " << std::endl;
  oclConnection :: Instance () -> print_ocl_objects ();
  std::cout << std::endl;
  
  std::cout << " * ocltest finished!" << std::endl;

	return true;

}
