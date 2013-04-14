# ifndef __OCL_FUNCTION_OBSERVER_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_FUNCTION_OBSERVER_HPP__
  
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclObservableDataObject.hpp"
  # include "oclSettings.hpp"
  
  
  
  
  /********************************
   ** class: oclFunctionObserver **
   ********************************/
  class oclFunctionObserver
  {
  

  
    public:
    
    
    
      /**
       * @name                  Constructors and destructors
       *                        Constructors and destructors.
       */
      //@{
      
      
      /**
       * @brief                 default constructor
       */
      oclFunctionObserver       (const int & max_num_args)
                               : mpp_args (
                                            (oclObservableDataObject ** const)
                                             malloc (max_num_args * sizeof (oclObservableDataObject *))
                                          ),
                                 m_num_args (0),
                                 m_max_num_args (max_num_args)
      {
      
        print_optional ("Ctor: oclFunctionObserver", v_level);
      
      }
      
      
      /**
       * @brief                 virtual destructor
       */
      virtual
      ~oclFunctionObserver      ( )
      {
      
        print_optional ("Dtor: oclFunctionObserver", v_level);
        
        /* free memory of array for registered objects */
        free (mpp_args);
      
      }
      
      
      //@}
      
      
      /**
       * @brief                 Register data objects add observer.
       *
       * @param  p_obj          Pointer to data object to be registered.
       */
      void
      register_obj              (oclObservableDataObject * const p_obj);
      
      
      /**
       * @brief                 Notify finished kernel.
       */
      void
      notify                    ( )
      const;
      

  
    private:
    
    
      /*****************
       ** member vars **
       *****************/
       
      /* registered objects */
            oclObservableDataObject ** const mpp_args;
      
      /* number of objects currently stored */
                                int          m_num_args;
      const                     int          m_max_num_args;
    
    
      /********************
       ** static members **
       ********************/
    
      /* private member for verbosity level of class */
      static const VerbosityLevel v_level;

  
  
  }; /* class oclFunctionObserver */
  
  
  
  
  /*************************************
   ** initialize static class members **
   *************************************/
  const VerbosityLevel oclFunctionObserver :: v_level = global_verbosity [OCL_FUNCTION_OBSERVER];
  
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  /**
   *                            Register data object.
   */
  void
  oclFunctionObserver ::
  register_obj                  (oclObservableDataObject * const p_obj)
  {
      
    print_optional ("oclObservableDataObject :: register_obj", v_level);
        
    /* check if there is space to store object */
    if (m_num_args < m_max_num_args)
    {
        
      /* add object, increment object counter */
      this -> mpp_args [m_num_args ++] = p_obj;
        
    }
    else
    {
        
      throw oclError ("Maximum number of stored objects reached!", "oclObservableDataObject :: register_obj");
        
    }
      
  }
  
  
  /**
   *                            Notify registered objects.
   */
  void
  oclFunctionObserver ::
  notify                        ( )
  const
  {
  
    print_optional ("oclObservableDataObject :: notify", v_level);
      
    /* loop over registered objects */
    for (int i = 0; i < m_num_args; i++)
    {
            
      /* call finish */
      mpp_args [i] -> finish ();
        
    }
    
  }
  
  
  
  
# endif /* __OCL_FUNCTION_OBSERVER_HPP__ */
