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
      oclGPUDataObject        (           T   * const   cpu_data,
                               const    int   &        num_elems)
                            : oclDataWrapper <T> (cpu_data, num_elems)
      {
      
        std::cout << "Ctor: \"oclGPUDataObject\"" << std::endl;
        
      }
      
      
      /**
       * @brief               "copy state" constructor
       */
      oclGPUDataObject        (                   T     * const  cpu_data,
                               const            int     &       num_elems,
                               const oclDataWrapper <T> &            obj)
                            : oclDataWrapper <T> (cpu_data, num_elems, obj)
      {
      
        std::cout << "Ctor: \"oclGPUDataObject\" ... copied state" << std::endl;
            
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
  
    std::cout << "oclGPUDataObject::prepare (" << oclDataObject :: getID () << ")" << std::endl;

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
    
//    std::cout << " -> keep data in GPU memory" << std::endl;

   // oclConnection :: Instance () -> getKernelArg (0, oclDataWrapper <T> :: mp_cpu_data, oclDataObject :: getSize ());
 
//    oclDataObject :: setGPUModified ();

//    oclDataWrapper <T> :: mp_cpu_data [0] = 1.5;
//    oclDataWrapper <T> :: mp_cpu_data [1] = 1.5;
//    oclDataWrapper <T> :: mp_cpu_data [2] = 1.5;
//    oclDataWrapper <T> :: mp_cpu_data [3] = 1.5;


    // update data state: available for use
    oclDataObject :: setUnlocked ();
  
//    std::cout << " /// getData !! ///" << std::endl;
//    getData ();
        
  }
  
  
  
  /**
   * @brief                   copy data to CPU memory
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
