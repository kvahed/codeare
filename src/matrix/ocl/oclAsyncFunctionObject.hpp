# ifndef __OCL_ASYNC_FUNCTION_OBJECT_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_ASYNC_FUNCTION_OBJECT_HPP__




  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclFunctionObserver.hpp"
  
  // pthread
  # include <pthread.h>
    
  
  
  
  /***********************************
   ** class: oclAsyncFunctionObject **
   **   (interface)                 **
   ***********************************/
  class oclAsyncFunctionObject
  {
  
  
  
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      bool                  m_active;
      oclFunctionObserver   m_observer;
      pthread_t             m_thread;
      
    
      
      /**********************
       ** member functions **
       **********************/
      
      /**
       * @name                  Constructors and destructors
       *                        Constructors and destructors.
       */
      //@{
      
      
      /**
       * @brief                 Constructor.
       *
       * @param  num_args       Number of arguments the observer
       *                        should allocate memory for.
       *
       */
      oclAsyncFunctionObject    ( const int & num_args )
                               : m_active   ( false ),
                                 m_observer ( num_args )
      {
      
        print_optional ("Ctor: oclAsyncFunctionObject", v_level);
      
      }
      
      
      /**
       * @brief                 Virtual destructor.
       */
      virtual
      ~oclAsyncFunctionObject    ( )
      {
      
        print_optional ("Dtor: oclAsyncFunctionObject", v_level);
      
      }
      
      
      //@}
      
  
      virtual
      void
      run_async           ( ) = 0;
      
      virtual
      void
      activate            () = 0;
            
      
    public:
    
    
      virtual
      bool
      isActive            ( )
      const = 0;
      
      virtual
      void
      join                ( )
      const = 0;


    private:
    
      /********************
       ** static members **
       ********************/
    
      /* private member for verbosity level of class */
      static const VerbosityLevel v_level;
    
  
  
  }; // class oclAsyncFunctionObject
  
  

  /*************************************
   ** initialize static class members **
   *************************************/
  const VerbosityLevel oclAsyncFunctionObject :: v_level = global_verbosity [OCL_ASYNC_FUNCTION_OBJECT];
    
  
  
  
# endif // __OCL_ASYNC_FUNCTION_OBJECT_HPP__
