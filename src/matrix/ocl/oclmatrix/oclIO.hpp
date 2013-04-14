# ifndef __OCL_IO_HPP__

  
  
  /************
   ** makros **
   ************/
  # define __OCL_IO_HPP__
  


  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclMatrix.hpp"
  
  // std C++
  # include <iostream>

  
  
  /**
   * @brief             output streaming operator for oclMatrix <T>
   */
  template <class T>
  std::ostream &
  operator<<            (std::ostream & os, const oclMatrix <T> & mat)
  {
    
    // assure data is available on CPU
    mat.getData ();
    
    // print oclMatrix <T> (use function for Matrix <T>)
    print (mat, os);
        
    // return output stream
    return os;
    
  }
  
  
/*  inline static std::ostream&  
  print (const Matrix<bool>& M, std::ostream &os) {
    
    for (size_t i = 0; i < M.Dim(0); i++) {
        for(size_t j = 0; j < M.Dim(1); j++)
			printf ("%zu ", M(i,j));
        printf("\n");
    }
    
    return os;
    
  }
*/  
  
# endif // __OCL_IO_HPP__
