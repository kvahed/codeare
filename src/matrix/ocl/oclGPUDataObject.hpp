# ifndef __OCL_GPU_DATA_OBJECT_HPP__



  /************
   ** makros **
   ************/
  #define __OCL_GPU_DATA_OBJECT_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclSettings.hpp"
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
      
        print_optional ("Ctor: \"oclGPUDataObject\" (id: %d)", oclDataObject :: getID (), v_level);

      }
      
      
      /**
       * @brief               "copy state" constructor
       */
      oclGPUDataObject        (                              T     * const  cpu_data,
                               const                       int     &       num_elems,
                                                oclDataWrapper <T> &             obj,
                                     oclDataObject :: CopyMode             copy_mode = oclDataObject :: NO_BUFFER)
                            : oclDataWrapper <T> (cpu_data, num_elems, obj, copy_mode)
      {
      
        print_optional ("Ctor: \"oclGPUDataObject\" ... copied state (new id: %d)", oclDataObject :: getID (), v_level);
        print_optional ("old id: %d", obj.getID (), v_level);
            
      }
      
      
      /**
       * @brief               virtual destructor
       */
      virtual
      ~oclGPUDataObject       ()
      {
      
        while (oclDataObject :: getLockState ())
          usleep (5);
      
        print_optional ("Dtor: \"oclGPUDataObject\"", VERB_HIGH);

      }

      
      //@}

    
      /**
       * @brief               inherited (oclDataObject)
       */
      virtual
      double
      prepare                 ();
      
      
      /**
       * @brief               inherited (oclObservableDataObject)
       */
      virtual
      double
      finish                  ();
      
      
      /**
       * @brief               inherited (oclDataWrapper)
       */
      virtual
      double
      getData                 ();



    private:
    
      /* private member for verbosity level of class */
      static const VerbosityLevel v_level;
    
    
     
  }; // class oclGPUDataObject
  


  /*************************************
   ** initialize static class members **
   *************************************/
  template <class T>
  const VerbosityLevel oclGPUDataObject <T> :: v_level = global_verbosity [OCL_GPU_DATA_OBJECT];

  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  
  /**
   * @brief                   prepare data object for use on GPU
   *                           -- load data to GPU if needed --
   *                          !! precondition: data not in use !!
   */
  template <class T>
  double
  oclGPUDataObject <T> ::
  prepare                     ()
  {

    print_optional ("oclGPUDataObject::prepare (%d)", oclDataObject :: getID (), v_level);
  
    double mem_time = .0;
    
    if (! oclDataObject :: getLockState ( ))
    {
  
      // set status: calculating (set available via finish ())
      oclDataObject :: setLocked ();
      
      // synchronize GPU data / load to GPU
      mem_time = oclDataWrapper <T> :: loadToGPU ();
    
      // notify modification of GPU data
      oclDataObject :: setGPUModified ();

    }

    return mem_time;
    
  }
  
  
  
  /**
   * @brief                   keep data in GPU memory (since it's a GPU object)
   */
  template <class T>
  double
  oclGPUDataObject <T> ::
  finish                      ()
  {

    print_optional ("oclGPUDataObject::finish (%d)", oclDataObject :: getID (), v_level);

    // time for memory transfer
    double mem_time = .0;
    
    /* check if buffer was grabbed */
    if (oclDataObject :: m_release_buffer)
    {
    
      print_optional (" oclGPUDataObject :: finish -> releaseBuffer", v_level);
      
      /* copy data to CPU */
      mem_time = oclDataWrapper <T> :: loadToCPU ( );
      
      /* release buffer */
      oclDataObject :: releaseBuffer ( );
    
    }

    // update CPU data
//    oclConnection :: Instance () -> loadToCPU (oclDataObject :: mp_gpu_buffer, oclDataWrapper <T> :: mp_cpu_data, oclDataObject :: m_size);

    // update data state: available for use
    oclDataObject :: setUnlocked ();
    
    return mem_time;
    
  }
  
  
  
  /**
   * @brief                   copy data to CPU memory
   */
  template <class T>
  double
  oclGPUDataObject <T> ::
  getData                     ()
  {
  
    print_optional ("oclGPUDataObject::getData", v_level);
    
    // memory transfer time
    double mem_time = .0;
    
    // check wether data is available or used on GPU
    if (oclDataObject :: getLockState ())
    {

      print_optional ("oclGPUDataObject::getData -> locked!", v_level);
      
      /* throw error */
      throw oclError ("Calculating on GPU ... data not available!", "oclGPUDataObject :: getData");
    
    }
    else
    {
    
      print_optional ("oclGPUDataObject::getData -> not locked", v_level);
    
      try
      {

        // synchronize CPU data with GPU
        mem_time = oclDataWrapper <T> :: loadToCPU ();

      }
      catch (const oclError & err)
      {
      
        print_optional (oclError (err, "oclGPUDataObject :: getData"), VERB_LOW);
      
      }

    }
    
    return mem_time;
    
  }
  
  
  
# endif // __OCL_GPU_DATA_OBJECT_HPP__
