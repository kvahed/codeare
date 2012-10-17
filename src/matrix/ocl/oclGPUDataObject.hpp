# ifndef __OCL_GPU_DATA_OBJECT_HPP__



  #define __OCL_GPU_DATA_OBJECT_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclDataWrapper.hpp"
  
  
  
  /*****************************
   ** class: oclGPUDataObject **
   **   (derived)             **
   *****************************/
  template <class T>
  class oclGPUDataObject : public oclDataWrapper <T>
  {
  
    
    public:
        
      /**
       * @name  constructors and destructors
       */
      //@{
      
      /**
       * @brief         constructor
       */
      oclGPUDataObject  (T * const cpu_data, const size_t & size)
                      : oclDataWrapper <T> (cpu_data, size)
      {
      
        std::cout << "Ctor: \"oclGPUDataObject\"" << std::endl;
        
        /* TODO */
        
      }
      
      //@}
    
      /**
       * @brief         inherited (oclDataObject)
       */
      virtual
      oclError
      prepare           () { /* TODO */ };
      
      /**
       * @brief         inherited (oclDataObject)
       */
      virtual
      oclError
      finish            () { /* TODO */ };
      
      /**
       * @brief         inherited (oclDataWrapper)
       */
      virtual
      T * const
      getData           () { /* TODO */ };
    
  
  }; // class oclGPUDataObject
  
  
  
# endif // __OCL_GPU_DATA_OBJECT_HPP__
