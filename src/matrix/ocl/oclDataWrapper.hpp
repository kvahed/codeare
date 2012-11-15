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
      oclDataWrapper        (T * const cpu_data, const size_t size)
                          : oclDataObject (size),
                            mp_cpu_data   (cpu_data)
      {
      
        std::cout << "Ctor: \"oclDataWrapper\"" << std::endl;
                        
      }
      
      
      /**
       * @brief             "copy state" constructor
       *
       */
      oclDataWrapper        (                   T     * const cpu_data,
                             const         size_t                 size,
                             const oclDataWrapper <T> &            obj)
                          : oclDataObject   (obj),
                            mp_cpu_data     (cpu_data)
      {
      
        std::cout << "Ctor: \"oclDataWrapper\" ... copied state" << std::endl;
              
        if (size != obj.getSize ())
        {
        
          throw " *!* Error: Sizes don't match! *!*";
        
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
        std::cout << "copy" << std::endl;
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
  
  
  
# endif // __OCL_DATA_WRAPPER_HPP__
