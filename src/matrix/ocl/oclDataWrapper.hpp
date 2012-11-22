# ifndef __OCL_DATA_WRAPPER_HPP__



  /************
   ** makros **
   ************/
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
       * @name              constructors and destructors
       */
      //@{
      
      
      /**
       * @brief             construtor
       *
       * @param             cpu_data    data in linear memory in RAM
       * @param             size        ... in bytes
       *
       */
      oclDataWrapper        (        T * const   cpu_data,
                             const int          num_elems)
                          : oclDataObject (num_elems, num_elems * sizeof (T)),
                            mp_cpu_data   (cpu_data)
      {
      
        std::cout << "Ctor: \"oclDataWrapper\"" << std::endl;
                        
      }
      
      
      /**
       * @brief             "copy state" constructor
       *
       */
      oclDataWrapper        (                   T     * const    cpu_data,
                             const            int               num_elems,
                                   oclDataWrapper <T> &               obj,
                                             bool             keep_buffer = false)
                          : oclDataObject   (obj, keep_buffer),
                            mp_cpu_data     (cpu_data)
      {
      
        std::cout << "Ctor: \"oclDataWrapper\" ... copied state" << std::endl;

        /* check if sizes match (copying allowed) */
        if (num_elems != obj.getNumElems())
        {
        
          throw " *!* Error: Num_elems don't match! *!*";
        
        }
        
        /* copy buffer of obj on GPU */
        if (keep_buffer && obj.bufferCopyable ())
        {
         
          // create buffer object
          oclConnection :: Instance () -> createBuffer (mp_cpu_data, oclDataObject :: getSize (), oclDataObject :: getID ());
    
          // update memory state
          oclDataObject :: setLoaded ();

        }
        else
        {
        
          /* data on GPU is newer than on CPU ? */
          if (obj.getMemState () && ! obj.getSyncState ())
          {
        
            // copy data of obj to CPU (so that data can be copied by oclMatrix)
            obj.loadToCPU ();

          }
        
        }
              
      }
      
      
      /**
       * @brief             virtual destructor
       */
      virtual
      ~oclDataWrapper       ()
      {
      
        std::cout << "Dtor: \"oclDataWrapper\"" << std::endl;
                
      }
        
      
      //@}
    
    
      /**
       * @brief             getter to cpu_data
       */
      virtual
      oclError &
      getData               () = 0;
      
      
      /**
       * @brief             print object's state to command line
       */
      virtual
      void
      print                 ();
      
      
      
    protected:

      
      /**
       * @brief             synchronize data on GPU with CPU
       *
       * @see               oclDataObject :: loadToGPU ()
       */
      virtual
      oclError &
      loadToGPU             ();
      
      
      /**
       * @brief             synchronize data on CPU with GPU
       *
       * @see               oclDataObject :: loadToCPU ()
       */
      virtual
      oclError &
      loadToCPU             ();
      
      
      /**********************
       ** member variables **
       **   (protected)    **
       **********************/
      T * const mp_cpu_data;            // pointer to cpu memory
      
      
    
  }; // class oclDataWrapper
  
  
  
  /**************************
   ** function definitions **
   **************************/


  /**
   * @brief                 -- refer to class definition --
   */
  template <class T>
  oclError &
  oclDataWrapper <T> ::
  loadToGPU                 ()
  {
      
    std::cout << "loadToGPU ()" << std::endl;

    // buffer on GPU exists ?
    if (! oclDataObject :: getMemState ())
    {

      // create buffer object
      oclConnection :: Instance () -> createBuffer (mp_cpu_data, oclDataObject :: getSize (), oclDataObject :: getID ());
    
      // update memory state
      oclDataObject :: setLoaded ();
    
    }

    if (oclDataObject :: mp_modified [CPU])
    {
      
      // update GPU data
      oclConnection :: Instance () -> loadToGPU (mp_cpu_data, oclDataObject :: getSize (), oclDataObject :: getBuffer ());

      // update states
      oclDataObject :: setSync ();

    }
    
    
  }
  
  
  
  /**
   * @brief                 -- refer to class definition --
   */
  template <class T>
  oclError &
  oclDataWrapper <T> ::
  loadToCPU                 ()
  {

    // buffer on GPU exists
    if (oclDataObject :: getMemState ())
    {
  
      std::cout << "loadToCPU ()" << std::endl;
          
      if (oclDataObject :: mp_modified [GPU])
      {
        
        // update CPU data
        oclConnection :: Instance () -> loadToCPU (oclDataObject :: mp_gpu_buffer, mp_cpu_data, oclDataObject :: m_size);
        
        // update states
        oclDataObject :: setSync ();
        
      }
      
    }
    else // ERROR if no buffer exists
    {
    
      std::cout << " *!* Error: no buffer on GPU (id:" << oclDataObject :: getID () << ") *!*" << std::endl;
    
    }
  
  }
  
  
  
  /**
   * @brief             print object's state to command line
   */
  template <class T>
  void
  oclDataWrapper <T> ::
  print                 ()
  {
  
    /* call method from super class */
    oclDataObject :: print ();
    
    /* add own state */
    std::cout << " *%* -oclDataWrapper-" << std::endl;
    std::cout << " *%*  -> data: ";
    if (oclDataObject :: getNumElems () > 10)
      std::cout << " << too large >> " << std::endl;
    else
    {
      for (int i = 0; i < oclDataObject :: getNumElems (); i++)
        std::cout << mp_cpu_data [i] << " ";
      std::cout << std::endl;
    }
  
  }

  
  
# endif // __OCL_DATA_WRAPPER_HPP__
