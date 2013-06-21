# ifndef __OCL_ASYNC_KERNEL_OBJECT_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_ASYNC_KERNEL_OBJECT_HPP__




  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclAsyncFunctionObject.hpp"  
  
  
  
  
  /*********************************
   ** class: oclAsyncKernelObject **
   **   (derived)                 **
   *********************************/
  class oclAsyncKernelObject : public oclKernelObject,
                               public oclAsyncFunctionObject
  {
  
  
    
    public:
    
    
      /**
       * @name                        Constructors and destructors
       *                              Constructors and destructors.
       */
      //@{
      
      
      /**
       * @brief                       Constructor (such as oclKernelObject).
       *
       * @param  kernel_name          Name of kernel to call.
       * @param      pp_args          Pointers to kernel's arguments.
       * @param     num_args          Number of kernel arguments.
       *
       */
      oclAsyncKernelObject            ( const   std::string &               kernel_name,
                                              oclDataObject * const * const     pp_args,
                                        const           int &                  num_args )
                                     : oclKernelObject        ( kernel_name, pp_args, num_args ),
                                       oclAsyncFunctionObject ( num_args )
      {
      
        print_optional ("Ctor: oclAsyncKernelObject", oclAsyncKernelObject :: v_level);
      
      }
      
      
      /**
       * @brief                       Virtual destructor.
       */
      virtual
      ~oclAsyncKernelObject           ( )
      {
      
        print_optional ("Dtor: oclAsyncKernelObject", oclAsyncKernelObject :: v_level);
      
      }
      
      
      //@}
    
    
      
      /**
       * @name                        Public getters.
       */
      //@{
      
      
      /**
       * @brief                       Get object's activity state.
       */
      inline
      virtual
      bool
      isActive                        ( )
      const;
      
      
      //@}
      
      
      
      /**
       * @brief                       Wait for asynchronous execution to finish.
       */
      virtual
      void
      join                            ( )
      const;
      
      
      
      /**
       * @brief                       Run kernel (inherited from oclKernelObject).
       */
      virtual
      void
      run                             ( );
    


    private:


      /**********************
       ** member functions **
       **   (private)      **
       **********************/

      
      /**
       * @brief                       Run kernel in asynchronous mode.
       */
      virtual
      void
      run_async                       ( );


      /**
       * @name                        Private setters.
       */
      //@{
      
      
      /**
       * @brief                       Set object active.
       */
      inline
      virtual
      void
      activate                        ( );
      
      
      //@}


      /**
       * @brief                       Function to be called by thread.
       */
      static
      void *
      call_thread                     ( void * obj )
      {
      
        /* cast to oclAsyncKernelObject */
        oclAsyncKernelObject * kernel_obj = (oclAsyncKernelObject *) obj;

        print_optional ("oclAsyncFunctionObject :: call_thread", v_level);

        /* run algorithm */
        kernel_obj -> run_async ( );

        return obj;
      
      }


      /********************
       ** static members **
       ********************/
    
      /* private member for verbosity level of class */
      static const VerbosityLevel v_level;


  
  }; /* oclAsyncKernelObject */
  


  /*************************************
   ** initialize static class members **
   *************************************/
  const VerbosityLevel oclAsyncKernelObject :: v_level = global_verbosity [OCL_ASYNC_KERNEL_OBJECT];
  
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  /**
   *                                  Get object's activity state.
   */
  inline
  bool
  oclAsyncKernelObject ::
  isActive                            ( )
  const
  {
  
    print_optional ("oclAsyncKernelObject :: isActive", oclAsyncKernelObject :: v_level);
    
    return oclAsyncFunctionObject :: m_active;
  
  }
  
  
  /**
   *                                  Set object active.
   */
  inline
  void
  oclAsyncKernelObject ::
  activate                            ( )
  {
  
    print_optional ("oclAsyncKernelObject :: activate", oclAsyncKernelObject :: v_level);
  
    /* set activity flag */
    oclAsyncFunctionObject :: m_active = true;
  
  }
  
  
  /**
   *                                  Run kernel in asynchronous mode.
   */
  void
  oclAsyncKernelObject ::
  run_async                           ( )
  {
  
    print_optional ("oclAsyncKernelObject :: run_async", oclAsyncKernelObject :: v_level);

    // run kernel
    cl::NDRange global_dims (256);
    cl::NDRange local_dims = cl::NullRange;
    oclConnection :: Instance () -> runKernel (global_dims, local_dims);

    /* deactivate asyncKernelObject */
    oclAsyncFunctionObject :: m_active = false;

    /* notify end of kernel execution */
    oclAsyncFunctionObject :: m_observer.notify ( );
  
  }

  
  
  /**
   *                                  Wait for async execution to finish.
   */
  void
  oclAsyncKernelObject ::
  join                                ( )
  const
  {
  
    print_optional ("oclAsyncKernelObject :: run", oclAsyncKernelObject :: v_level);
    
    /* call join on thread handle */
    pthread_join (oclAsyncFunctionObject :: m_thread, NULL);
  
  }
  
  
  
  /**
   *                                  Run kernel (inherited from oclKernelObject).
   */
  void
  oclAsyncKernelObject ::
  run                                 ( )
  {
  
    print_optional ("oclAsyncKernelObject :: run", oclAsyncKernelObject :: v_level);
    
    oclConnection * const oclCon = oclConnection :: Instance ( );
    
    /* activate asyncKernelObject */
    this -> activate ( );
    
    // activate kernel
    oclCon -> activateKernel (m_kernel_name);
    
    // prepare kernel arguments (load to gpu)
    for (int i = 0; i < m_num_args; i++)
    {

      // prepare argument
      mpp_args [i] -> prepare ();

      // register argument at kernel
      oclCon -> setKernelArg (i, mpp_args [i]);

      try {

        /* register argument at function observer */
        oclAsyncFunctionObject :: m_observer.register_obj (mpp_args [i]);

      }
      catch (oclError & err)
      {
      
        throw oclError (err, "oclAsyncKernelObject :: run");
      
      }

    }

    /* start thread for asynchronous kernel call */
    pthread_create (& (oclAsyncFunctionObject :: m_thread), NULL, &oclAsyncKernelObject :: call_thread, this);

//    this -> join ( );

//    pthread_join (oclAsyncFunctionObject :: m_thread, NULL);

  }  
  
  
  
# endif /* __OCL_ASYNC_KERNEL_OBJECT_HPP__ */
