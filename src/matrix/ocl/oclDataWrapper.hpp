# ifndef __OCL_DATA_WRAPPER_HPP__



  # define __OCL_DATA_WRAPPER_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclDataObject.hpp"
  
  
  
  /***************************
   ** class: oclDataWrapper **
   **   (derived)           **
   ***************************/
  template <class T>
  class oclDataWrapper : public oclDataObject
  {
    
    
    public:
    
      /**
       * @brief         getter to cpu_data
       */
      virtual
      T * const
      getData           () const = 0;
    
    
    protected:
        
      // constructor
      oclDataWrapper    (T * const cpu_data, const size_t size)
                      : oclDataObject (),
                        mp_cpu_data   (cpu_data),
                        m_size        (size)
      {
      
        std::cout << "Ctor: \"oclDataWrapper\"" << std::endl;
        
        /* TODO */
        
      }
      
      // destructor
      virtual
      ~oclDataWrapper   ()
      {
        std::cout << "Dtor: \"oclDataWrapper\"" << std::endl;
        /* TODO */
      }
        
      // pointer to cpu memory
      T * const mp_cpu_data;
      
      // size of cpu_data
      size_t m_size;
    
    
  }; // class oclDataWrapper
  
  
  
# endif // __OCL_DATA_WRAPPER_HPP__
