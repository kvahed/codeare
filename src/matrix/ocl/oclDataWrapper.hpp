# ifndef __OCL_DATA_WRAPPER_HPP__



  /************
   ** makros **
   ************/
  # define __OCL_DATA_WRAPPER_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclSettings.hpp"
  # include "oclDataObject.hpp"
  
  // std C++ libraries
  # include <sstream>
  
  
  
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
      
        print_optional ("Ctor: \"oclDataWrapper\"", VERB_HIGH);
                        
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
      
        print_optional ("Ctor: \"oclDataWrapper\" ... copied state", VERB_HIGH);

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
      
        print_optional ("Dtor: \"oclDataWrapper\"", VERB_HIGH);
                
      }
        
      
      //@}
    
    
      /**
       * @brief             getter to cpu_data
       */
      virtual
      void
      getData               ()
      throw (oclError) = 0;
      
      
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
      void
      loadToGPU             ();
      
      
      /**
       * @brief             synchronize data on CPU with GPU
       *
       * @see               oclDataObject :: loadToCPU ()
       */
      virtual
      void
      loadToCPU             ()
      throw (oclError);
      
      
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
  void
  oclDataWrapper <T> ::
  loadToGPU                 ()
  {
      
    print_optional ("oclDataWrapper :: loadToGPU ()", VERB_HIGH);

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
  void
  oclDataWrapper <T> ::
  loadToCPU                 ()
  throw (oclError)
  {

    print_optional ("oclDataWrapper :: loadToCPU ()", VERB_HIGH);

    // buffer on GPU exists
    if (oclDataObject :: getMemState ())
    {
  
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
    
      stringstream tmp;
      tmp << "No buffer on GPU (id:" << oclDataObject :: getID () << ")";
      
      throw oclError (tmp.str (), "oclDataWrapper :: loadToCPU", oclError :: WARNING);
    
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
    if (oclDataObject :: getNumElems () > 16)
      std::cout << " << too large >> " << std::endl;
    else
    {
      for (int i = 0; i < oclDataObject :: getNumElems (); i++)
        std::cout << mp_cpu_data [i] << " ";
      std::cout << std::endl;
    }
  
  }

  
  
# endif // __OCL_DATA_WRAPPER_HPP__
