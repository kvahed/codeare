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
      
      /**
       * @brief         destructor
       */
      virtual
      ~oclGPUDataObject ()
      {
        std::cout << "Dtor: \"oclGPUDataObject\"" << std::endl;
        /* TODO */
      }
      
      //@}
    
      /**
       * @brief         inherited (oclDataObject)
       */
      virtual
      oclError &
      prepare           (const int num);
      
      /**
       * @brief         inherited (oclDataObject)
       */
      virtual
      oclError &
      finish            (const int num);
      
      /**
       * @brief         inherited (oclDataWrapper)
       */
      virtual
      T * const
      getData           ()
      const;

    
  
  }; // class oclGPUDataObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  template <class T>
  oclError &
  oclGPUDataObject <T> ::
  prepare               (const int num)
  {
  
    std::cout << "oclGPUDataObject::prepare" << std::endl;

    oclConnection :: Instance () -> setKernelArg (num, oclDataWrapper<T> :: mp_cpu_data, oclDataWrapper<T> :: m_size * sizeof (T), CL_MEM_READ_WRITE);

  }
  
  
  template <class T>
  oclError &
  oclGPUDataObject <T> ::
  finish                (const int num)
  {

    std::cout << "oclGPUDataObject::finish" << std::endl;
    
    oclConnection :: Instance () -> getKernelArg (num, oclDataWrapper<T> :: mp_cpu_data, oclDataWrapper<T> :: m_size * sizeof (T));
  
  }
  
  
  
  template <class T>
  T * const
  oclGPUDataObject <T> ::
  getData               ()
  const
  {
    /* TODO */
  }
  
  
  
# endif // __OCL_GPU_DATA_OBJECT_HPP__
