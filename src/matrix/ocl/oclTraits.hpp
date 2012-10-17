# ifndef __OCL_TRAITS_HPP__



  # define __OCL_TRAITS_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclDataWrapper.hpp"
  # include "oclGPUDataObject.hpp"
  
  
  
  /***********************
   ** struct: oclTraits **
   **   (base struct)   **
   ***********************/
  template <class T>
  struct oclTraits
  {
    /* -- */
  }; // struct oclTraits
  
  
  
  /*******************************
   ** struct: oclTraits <float> **
   **   (single precision)      **
   *******************************/
  template <>
  struct oclTraits <float>
  {
    
    
    /**********************
     ** type definitions **
     **********************/
    typedef float elem_type;
    
    
    /**
     * @brief         Create oclDataObject
     */
    static inline
    oclDataWrapper <elem_type> *
    make_GPU_Obj      (elem_type * const cpu_arg, const size_t & size)
    {
    
      std::cout << "make_GPU_Obj" << std::endl;
    
      return new oclGPUDataObject <elem_type> (cpu_arg, size);
      
    }
    
    
  };
  
  
  
  
# endif // __OCL_TRAITS_HPP__
