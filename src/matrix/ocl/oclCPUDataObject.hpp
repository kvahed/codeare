# ifndef __OCL_CPU_DATA_OBJECT_HPP__



  /************
   ** makros **
   ************/
  #define __OCL_CPU_DATA_OBJECT_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclSettings.hpp"
  # include "oclDataWrapper.hpp"
  
  
  
  /*****************************
   ** class: oclCPUDataObject **
   **   (derived)             **
   *****************************/
  template <class T>
  class oclCPUDataObject : public oclDataWrapper <T>
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
      oclCPUDataObject        (           T   * const   cpu_data,
                               const    int   &        num_elems)
                            : oclDataWrapper <T> (cpu_data, num_elems)
      {
      
        print_optional ("Ctor: \"oclCPUDataObject\"", VERB_HIGH);
        
      }
      
      
      /**
       * @brief               "copy state" constructor
       */
/*      oclCPUDataObject        (                   T     * const    cpu_data,
                               const            int     &         num_elems,
                                     oclDataWrapper <T> &               obj,
                                               bool             keep_buffer = false)
                            : oclDataWrapper <T> (cpu_data, num_elems, obj, keep_buffer)
      {
      
        print_optional ("Ctor: \"oclCPUDataObject\" ... copied state", VERB_HIGH);
            
      }
*/      
      
      /**
       * @brief               virtual destructor
       */
      virtual
      ~oclCPUDataObject       ()
      {
      
        print_optional ("Dtor: \"oclCPUDataObject\"", VERB_HIGH);

      }

      
      //@}

    
      /**
       * @brief               inherited (oclDataObject)
       */
      virtual
      void
      prepare                 ();
      
      
      /**
       * @brief               inherited (oclObservableDataObject)
       */
      virtual
      void
      finish                  ();
      
      
      /**
       * @brief               inherited (oclDataWrapper)
       */
      virtual
      void
      getData                 ()
      throw (oclError);

    
     
  }; // class oclCPUDataObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  
  /**
   * @brief                   prepare data object for use on GPU
   *                           -- load data to GPU if needed --
   *                          !! precondition: data not in use !!
   */
  template <class T>
  void
  oclCPUDataObject <T> ::
  prepare                     ()
  {

    print_optional ("oclCPUDataObject::prepare (%d)", oclDataObject :: getID (), VERB_MIDDLE);
  
    // set status: calculating (set available via finish ())
    oclDataObject :: setLocked ();
    
    // synchronize GPU data / load to GPU
    oclDataWrapper <T> :: loadToGPU ();
    
    // notify modification of GPU data
    oclDataObject :: setGPUModified ();

  }
  
  
  
  /**
   * @brief                   load data to CPU memory (since it's a CPU object)
   */
  template <class T>
  void
  oclCPUDataObject <T> ::
  finish                      ()
  {

    print_optional ("oclCPUDataObject::finish", VERB_HIGH);

    // update data state: available for use
    oclDataObject :: setUnlocked ();

    // copy data to CPU memory
    this -> getData ();

  }
  
  
  
  /**
   * @brief                   copy data to CPU memory
   */
  template <class T>
  void
  oclCPUDataObject <T> ::
  getData                     ()
  throw (oclError)
  {
  
    print_optional ("oclCPUDataObject::getData", VERB_HIGH);
    
    // check wether data is available or used on GPU
    if (oclDataObject :: getLockState ())
    {
    
      /* throw error */
      throw oclError ("Calculating on GPU ... data not available!", "oclCPU_DataObject :: getData");
    
    }
    else
    {
    
      try
      {

        // synchronize CPU data with GPU
        oclDataWrapper <T> :: loadToCPU ();

      }
      catch (const oclError & err)
      {
      
        print_optional (oclError (err, "oclCPUDataObject :: getData"), VERB_LOW);
      
      }

    }    
    
  }
  
  
  
# endif // __OCL_CPU_DATA_OBJECT_HPP__
