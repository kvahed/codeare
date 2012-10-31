# ifndef __OCL_GPU_DATA_OBJECT_HPP__



  /************
   ** makros **
   ************/
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
       * @name                constructors and destructors
       */
      //@{
      
      
      /**
       * @brief               constructor
       *
       * @param               cpu_data  @see oclDataWrapper <T>
       * @param               size      @see ocldataWrapper <T>
       *
       */
      oclGPUDataObject        (           T   * const  cpu_data,
                               const size_t   &        size)
                            : oclDataWrapper <T> (cpu_data, size)
      {
      
        std::cout << "Ctor: \"oclGPUDataObject\"" << std::endl;
        
      }
      
      
      /**
       * @brief               virtual destructor
       */
      virtual
      ~oclGPUDataObject       ()
      {
      
        std::cout << "Dtor: \"oclGPUDataObject\"" << std::endl;

      }

      
      //@}

    
      /**
       * @brief               inherited (oclDataObject)
       */
      virtual
      oclError &
      prepare                 ();
      
      
      /**
       * @brief               inherited (oclObservableDataObject)
       */
      virtual
      oclError &
      finish                  ();
      
      
      /**
       * @brief               inherited (oclDataWrapper)
       */
      virtual
      oclError &
      getData                 ();

    
  
  }; // class oclGPUDataObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  
  /**
   * @brief                   prepare data object for use on GPU
   *                           -- load data to GPU if needed --
   *                          !! precondition: data not in use !!
   */
  template <class T>
  oclError &
  oclGPUDataObject <T> ::
  prepare                     ()
  {
  
    // set status: calculating (set available via finish ())
    oclDataObject :: setLocked ();
  
    std::cout << "oclGPUDataObject::prepare" << std::endl;

    // synchronize GPU data / load to GPU
    oclDataWrapper <T> :: loadToGPU ();
    
    // notify modification of GPU data
    oclDataObject :: setGPUModified ();

  }
  
  
  
  /**
   * @brief                   keep data in GPU memory (since it's a GPU object)
   */
  template <class T>
  oclError &
  oclGPUDataObject <T> ::
  finish                      ()
  {

    std::cout << "oclGPUDataObject::finish" << std::endl;
    
    std::cout << " -> keep data in GPU memory" << std::endl;
    
    // update data state: available for use
    oclDataObject :: setUnlocked ();
        
  }
  
  
  
  /**
   * @brief                   copy data to CPU memory, return pointer
   */
  template <class T>
  oclError &
  oclGPUDataObject <T> ::
  getData                     ()
  {
  
    std::cout << "oclGPUDataObject::getData" << std::endl;
    
    // check wether data is available or used on GPU
    if (oclDataObject :: getLockState ())
    {
    
      std::cout << " *!* calculating on GPU ... data not available *!* " << std::endl;
    
    }
    else
    {
    
      // synchronize CPU data with GPU
      oclDataWrapper <T> :: loadToCPU ();

    }    
    
  }
  
  
  
# endif // __OCL_GPU_DATA_OBJECT_HPP__
