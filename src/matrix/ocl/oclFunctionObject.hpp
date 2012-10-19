# ifndef __OCL_FUNCTION_OBJECT_HPP__



  # define __OCL_FUNCTION_OBJECT_HPP__



  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclDataObject.hpp"
  
  
  
  /**************************
   ** forward declarations **
   **************************/
//  class oclDataObject;
  

  
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
      oclDataObject     * const * const mpp_args;
      int                               m_num_args;
  
      // constructor
      oclFunctionObject   (oclDataObject * const * const pp_args,
                           int                           num_args )
                        : mpp_args     (pp_args),
                          m_num_args (num_args)
      {
        std::cout << "Ctor: \"oclFunctionObject\"" << std::endl;
        /* TODO */
      }
      
      virtual
      ~oclFunctionObject ()
      {
        std::cout << "Dtor: \"oclFunctionObject\"" << std::endl;
        /* TODO */
      }
  
  
  }; // class oclFunctionObject
  
  
  
# endif // __OCL_FUNCTION_OBJECT_HPP__
