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
# include "ocl/test_ops.hpp"



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
oclInit2D              ( size_t dimX, size_t dimY, const int & offset = 0 )
{

  oclMatrix <T> ocl_mat (dimX, dimY);
  for (size_t i = 0; i < dimX; i++)
    for (size_t j = 0; j < dimY; j++)
      ocl_mat (i, j) = (i + j + offset) / ((dimX + dimY) / 2);
      
  return ocl_mat;

}


template <class T>
Matrix <T>
Init2D              ( size_t dimX, size_t dimY, const int & offset = 0 )
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


template <class T, class S>
bool
oclComparisonsTestC   ( bool verbose )
{

  std::cout <<                          std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout << " ** comparisons ** " << std::endl;
  std::cout << " ***************** " << std::endl;
  std::cout <<                          std::endl;

  size_t  dimX     = 2048,
          dimY     = 2048;
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl << std::endl;

  /* compare matrix and scalar */
  {
    std::cout << "\t\t  * Matrix: A, scalar: s * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
                  S      scalar = (S) 3;
             Matrix <T>     mat = Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    // Tests //
    if (verbose) std::cout << " * A == s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opEqual <T, S, bool> ());
    if (verbose) std::cout << " * A != s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opInequal <T, S, bool> ());
/*    if (verbose) std::cout << " * s == A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opEqual <S, T, bool> ());
    if (verbose) std::cout << " * s != A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opInequal <S, T, bool> ());
*/  }
  
  /* compare matrices */
  {
    std::cout <<                                      std::endl;
    std::cout << "\t\t  * Matrix: A, Matrix: B * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
    const       int pos_x = 10,
                    pos_y = 11;
             Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 = Init2D <S> (dimX, dimY);
                            mat1 (pos_x, pos_y) += (S) 1;
          oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S>    ocl_mat2 = oclInit2D <S> (dimX, dimY);
                           ocl_mat1 (pos_x, pos_y) += (S) 1;
    // Tests //
    if (verbose) std::cout << " * A == B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opEqual <T, S, bool> ());
    if (verbose) std::cout << " * A != B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opInequal <T, S, bool> ());
  }

}


template <class T, class S>
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
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl << std::endl;

  /* compare matrix and scalar */
  {
    std::cout << "\t\t  * Matrix: A, scalar: s * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
                  S      scalar = (S) 3;
             Matrix <T>     mat = Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    // Tests //
    if (verbose) std::cout << " * A == s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opEqual <T, S, bool> ());
    if (verbose) std::cout << " * A != s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opInequal <T, S, bool> ());
    if (verbose) std::cout << " * A >  s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opGreater <T, S, bool> ());
    if (verbose) std::cout << " * A >= s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opGEqual <T, S, bool> ());
    if (verbose) std::cout << " * A <  s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opLess <T, S, bool> ());
    if (verbose) std::cout << " * A <= s                      ";
    operator_eval_MS <T, S, bool> (mat, ocl_mat, scalar, opLEqual <T, S, bool> ());
    if (verbose) std::cout << " * s == A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opEqual <S, T, bool> ());
    if (verbose) std::cout << " * s != A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opInequal <S, T, bool> ());
    if (verbose) std::cout << " * s >  A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opGreater <S, T, bool> ());
    if (verbose) std::cout << " * s >= A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opGEqual <S, T, bool> ());
    if (verbose) std::cout << " * s <  A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opLess <S, T, bool> ());
    if (verbose) std::cout << " * s <= A                      ";
    operator_eval_SM <S, T, bool> (scalar, mat, ocl_mat, opLEqual <S, T, bool> ());
  }
  
  /* compare matrices */
  {
    std::cout <<                                      std::endl;
    std::cout << "\t\t  * Matrix: A, Matrix: B * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
    const       int pos_x = 10,
                    pos_y = 11;
             Matrix <T>     mat1 = Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 = Init2D <S> (dimX, dimY);
                            mat1 (pos_x, pos_y) += (S) 1;
          oclMatrix <T>    ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S>    ocl_mat2 = oclInit2D <S> (dimX, dimY);
                           ocl_mat1 (pos_x, pos_y) += (S) 1;
    // Tests //
    if (verbose) std::cout << " * A == B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opEqual <T, S, bool> ());
    if (verbose) std::cout << " * A != B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opInequal <T, S, bool> ());
    if (verbose) std::cout << " * A >  B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opGreater <T, S, bool> ());
    if (verbose) std::cout << " * A >= B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opGEqual <T, S, bool> ());
    if (verbose) std::cout << " * A <  B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opLess <T, S, bool> ());
    if (verbose) std::cout << " * A <= B                      ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opLEqual <T, S, bool> ());
    if (verbose) std::cout << " * C = A && B                  ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opAND <T, S, bool> ());
    if (verbose) std::cout << " * C = A || B                  ";
    operator_eval_MM <T, S, bool> (mat1, mat2, ocl_mat1, ocl_mat2, opOR <T, S, bool> ());
  }

}


template <class T, class S>
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
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl << std::endl;

  /* M = op (M, s);
     M = op (s, M); */
  {
    std::cout << "\t\t  * Matrix: A, scalar: s * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
                  S      scalar = (S) 33.3;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    // Tests //
    if (verbose) std::cout << " * M = A + s                   ";
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opAdd <T, S, T> ());
    if (verbose) std::cout << " * M = A - s                   ";
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opSub <T, S, T> ());
    if (verbose) std::cout << " * M = A * s                   ";
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opMult <T, S, T> ());
    if (verbose) std::cout << " * M = A / s                   ";
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opDiv <T, S, T> ());
  }

  /* decrement matrix (uniform) */
  {
    if (verbose) std::cout << " * M -= s                      ";
    const         S      scalar = (S) 33.3;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opSubAss <T, S, T> ());
  }
  
  /* increment matrix (uniform) */
  {
    if (verbose) std::cout << " * M += s                      ";
    const         S      scalar = (S) 33.3;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opAddAss <T, S, T> ());
  }

  /* elementwise multiplacation of matrix and scalar */
  {
    if (verbose) std::cout << " * C *= s                      ";
    const         S      scalar = (S) (10 * (1./3.));
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat, ocl_mat, scalar, opMultAss <T, S, T> ());   
  }
 
  /* M = op (M1, M2); */
  {
    std::cout <<                                      std::endl;
    std::cout << "\t\t  * Matrix: A, Matrix: B * " << std::endl;
    std::cout <<                                      std::endl;
       Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
       Matrix <S>     mat2 =    Init2D <S> (dimX, dimY);
    oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, dimY);
    if (verbose) std::cout << " * M = A + B                   ";
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opAdd <T, S, T> ());
    if (verbose) std::cout << " * M = A - B                   ";
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opSub <T, S, T> ());
    if (verbose) std::cout << " * M = A .* B                  ";
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opMult <T, S, T> ());
    if (verbose) std::cout << " * M = A ./ B                  ";
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opDiv <T, S, T> ());
  }

  /* increment matrix (non uniform) */
  {
    if (verbose) std::cout << " * A += B                      ";
             Matrix <T>         mat =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat_inc =    Init2D <S> (dimX, dimY);
          oclMatrix <T>     ocl_mat = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat_inc = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat, mat_inc, ocl_mat, ocl_mat_inc, opAddAss <T, S, T> ());
  }

  /* decrement matrix (non uniform) */
  {
    if (verbose) std::cout << " * A -= B                      ";
             Matrix <T>         mat =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat_dec =    Init2D <S> (dimX, dimY);
          oclMatrix <T>     ocl_mat = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat_dec = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat, mat_dec, ocl_mat, ocl_mat_dec, opSubAss <T, S, T> ());
  }

  /* elementwise multiplacation of two matrices */
  {
    if (verbose) std::cout << " * C .*= A                     ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 =    Init2D <S> (dimX, dimY);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat2  = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opMultAss <T, S, T> ());   
  }

  /* elementwise division of two matrices */
  {
    if (verbose) std::cout << " * C ./= A                     ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 =    Init2D <S> (dimX, dimY);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opDivAss <T, S, T> ());
  }

}



template <class T, class S>
bool
oclArithmeticsTestExt   ( bool verbose )
{

  std::cout <<                                std::endl;
  std::cout << " *********************** " << std::endl;
  std::cout << " ** arithmetics (ext) ** " << std::endl;
  std::cout << " *********************** " << std::endl;
  std::cout <<                                std::endl;

  size_t  dimX     = 2048,
          dimY     = 2048;
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl << std::endl;

  /* M = op (M, s);
     M = op (s, M); */
  {
    std::cout << "\t\t  * Matrix: A, scalar: s * " << std::endl;
    std::cout <<                                      std::endl;
    // Init //
                  S      scalar = (S) 33.3;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    // Tests //
    if (verbose) std::cout << " * M = s / A                   ";
    operator_eval_SM <S, T, T> (scalar, mat, ocl_mat, opDiv <S, T, T> ());
    if (verbose) std::cout << " * M = s + A                   ";
    operator_eval_SM <S, T, T> (scalar, mat, ocl_mat, opAdd <S, T, T> ());
    if (verbose) std::cout << " * M = s - A                   ";
    operator_eval_SM <S, T, T> (scalar, mat, ocl_mat, opSub <S, T, T> ());
    if (verbose) std::cout << " * M = s * A                   ";
    operator_eval_SM <S, T, T> (scalar, mat, ocl_mat, opMult <S, T, T> ());
  }

  /* elementwise division of matrix by scalar */
  {
    if (verbose) std::cout << " * C /= s                      ";
    const         S       scalar = (S) 3.33;
             Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat1, ocl_mat1, scalar, opDivAss <T, S, T> ());
  }

  /* raise to higher power */
  {
    if (verbose) std::cout << " * C = A ^ p                   ";
    const         S           p = (S) 1.55;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat, ocl_mat, p, opPow <T, S, T> ());
  }

  /* raise to higher power */
  {
    if (verbose) std::cout << " * M ^= p                      ";
    const         S           p = (S) 1.55;
             Matrix <T>     mat =    Init2D <T> (dimX, dimY);
          oclMatrix <T> ocl_mat = oclInit2D <T> (dimX, dimY);
    operator_eval_MS <T, S, T> (mat, ocl_mat, p, opPowAss <T, S, T> ());
  }
  
  /* unary minus */
  {
    if (verbose) std::cout << " * C = -A                      ";
            Matrix <T>      mat1                 (dimX, dimY);
    const   Matrix <S>      mat2 =    Init2D <S> (dimX, dimY);
          oclMatrix <T> ocl_mat1                 (dimX, dimY);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opMinus <T, S, T> ());
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




template <class T, class S>
bool
oclLinalgTest           ( bool verbose )
{

  std::cout <<                             std::endl;
  std::cout << " ******************** " << std::endl;
  std::cout << " ** Linear Algebra ** " << std::endl;
  std::cout << " ******************** " << std::endl;
  std::cout <<                             std::endl;

  size_t  dimX     = 2048,
          dimY     = 1024;
  
  std::cout << " size: " << dimX << " x " << dimY << std::endl << std::endl;

  /* transpose matrix */
  {
    if (verbose) std::cout << " * C = ! A                     ";
             Matrix <T>     mat1                 (dimX, dimY);
    const    Matrix <S>     mat2 =    Init2D <S> (dimX, dimY);
          oclMatrix <T> ocl_mat1                 (dimX, dimY);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, dimY);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opTrans <T, S, T> ()); 
  }

  /* matrix product */
  {
    if (verbose) std::cout << " * C = A * B                   ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 =    Init2D <S> (dimY, dimX);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimY, dimX);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opMatProd <T, S, T> ());
  }
  
  /* matrix vector product */
  {
    if (verbose) std::cout << " * y = A * x                   ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, dimY);
    const    Matrix <S>     mat2 =    Init2D <S> (dimY,    1);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, dimY);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimY,    1);
    operator_eval_MM <T, S, T> (mat1, mat2, ocl_mat1, ocl_mat2, opMatProd <T, S, T> ());
  }

}



template <class T, class S>
bool
oclLinalgTestExt            ( bool verbose )
{

  size_t  dimX     = 2048;

  double time_s, time_ocl;
  MyTimer t;
  
  std::cout << std::endl << " size: " << dimX << " x 1" << std::endl << std::endl;
  
  /* dot product */
  {
    if (verbose) std::cout << " * d = <u, v>                 ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, 1);
    const    Matrix <S>     mat2 =    Init2D <S> (dimX, 1);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, 1);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, 1);
               t.tic (time_default);
    T dp_s   =     mat1.dot (    mat2);
    time_s   = t.tic (time_default);
    T dp_ocl = ocl_mat1.dot (ocl_mat2);
    time_ocl = t.tic (time_default);
    assert (dp_s == dp_ocl);
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

}



template <class T, class S>
bool
oclLinalgTestExtC          ( bool verbose )
{

  size_t  dimX     = 2048;
  
  double time_s, time_ocl;
  MyTimer t;

  std::cout << std::endl << " size: " << dimX << " x 1" << std::endl << std::endl;
  
  /* complex dot product */
  {
    if (verbose) std::cout << " * d = <u, v>                 ";
             Matrix <T>     mat1 =    Init2D <T> (dimX, 1);
    const    Matrix <S>     mat2 =    Init2D <S> (dimX, 1);
          oclMatrix <T> ocl_mat1 = oclInit2D <T> (dimX, 1);
    const oclMatrix <S> ocl_mat2 = oclInit2D <S> (dimX, 1);
               t.tic (time_default);
    T dp_s   =     mat1.dotc (    mat2);
    time_s   = t.tic (time_default);
    T dp_ocl = ocl_mat1.dotc (ocl_mat2);
    time_ocl = t.tic (time_default);
    assert (dp_s == dp_ocl);
    if (verbose) std::cout << "passed! (time_s: " << time_s << "s, time_ocl: " << time_ocl << "s)" << std::endl;
  }

}



template <class T, class S>
void
exec_tests      (bool verbose)
{

  oclCreatorsTest       <T>    (verbose);
  oclArithmeticsTest    <T, S> (verbose);
  oclArithmeticsTestExt <T, S> (verbose);
  oclLinalgTest         <T, S> (verbose);
  oclLinalgTestExt      <T, S> (verbose);
  oclAccessTest         <T>    (verbose);
  oclComparisonsTest    <T, S> (verbose);

}



template <class T, class S>
void
exec_testsC     (bool verbose)
{

  oclCreatorsTest     <T>    (verbose);
  oclArithmeticsTest  <T, S> (verbose);
  oclLinalgTest       <T, S> (verbose);
  oclLinalgTestExtC   <T, S> (verbose);
  oclAccessTest       <T>    (verbose);
  oclComparisonsTestC <T, S> (verbose);

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
/*    std::cout << std::endl << std::endl;
    std::cout << " ----------- // float, float  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <float, float>  (true);
    std::cout << std::endl << std::endl;
*//*    std::cout << " ----------- // float, double  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <float, double>  (true);
    std::cout << std::endl << std::endl;
*//*    std::cout << " ----------- // double, double \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <double, double> (true);
    std::cout << std::endl << std::endl;
*//*    std::cout << " ----------- // double, float  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_tests <double, float>  (true);
    std::cout << std::endl;
*/
    std::cout << std::endl << std::endl;
    std::cout << " ----------- // cxfl, cxfl  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_testsC <cxfl, cxfl> (true);
    std::cout << std::endl << std::endl;
    std::cout << " ----------- // cxdb, cxdb  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_testsC <cxdb, cxdb> (true);
    std::cout << std::endl << std::endl;
/*    std::cout << std::endl << std::endl;
    std::cout << " ----------- // cxfl, float  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_testsC <cxfl, float> (true);
    std::cout << std::endl << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << " ----------- // cxdb, double  \\\\ ------------ " << std::endl;
    std::cout << std::endl;
    exec_testsC <cxdb, double> (true);
    std::cout << std::endl << std::endl;
*/  }
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
