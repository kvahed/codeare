# ifndef __OCL_FUNCTION_OBJECT_HPP__



  # define __OCL_FUNCTION_OBJECT_HPP__



  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclDataObject.hpp"
  

  
  /******************************
   ** class: oclFunctionObject **
   **   (abstract)             **
   ******************************/
  class oclFunctionObject
  {
  
  
    public:
    
      virtual
      void
      run () = 0;
  
  
    protected:
    
      // function arguments
      oclDataObject * args;
  
  
  }; // class oclFunctionObject
  
  
  
# endif // __OCL_FUNCTION_OBJECT_HPP__
