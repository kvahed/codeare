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
       * @param         args          Pointers to kernel's arguments.
       * @param     num_args          Number of kernel arguments.
       *
       */
      oclAsyncKernelObject            ( const   std::string &               kernel_name,
                                              oclDataObject * const * const     pp_args,
                                        const           int &                  num_args )
                                     : oclKernelObject        ( kernel_name, pp_args, num_args ),
                                       oclAsyncFunctionObject ( num_args )
      {
      
        print_optional ("CTOR: oclAsyncKernelObject", oclAsyncKernelObject :: v_level);
      
      }
      
      
      /**
       * @brief                       Virtual destructor.
       */
      virtual
      ~oclAsyncKernelObject           ( )
      {
      
        print_optional ("DTOR: oclAsyncKernelObject", oclAsyncKernelObject :: v_level);
      
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
       * @brief                       Run kernel (inherited from oclKernelObject).
       */
      virtual
      void
      run                             ( )
      const;
    


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
      run_async                       ( )
      const;


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
  const
  {
  
    print_optional ("oclAsyncKernelObject :: run_async", oclAsyncKernelObject :: v_level);
  
    /* TODO */
  
  }
  
  
  /**
   *                                  Run kernel (inherited from oclKernelObject).
   */
  void
  oclAsyncKernelObject ::
  run                                 ( )
  const
  {
  
    print_optional ("oclAsyncKernelObject :: run", oclAsyncKernelObject :: v_level);
    
    /* TODO */
    
  }
  
  
  
  
# endif /* __OCL_ASYNC_KERNEL_OBJECT_HPP__ */
