# ifndef __OCL_ASYNC_FUNCTION_OBJECT_HPP__



  /************
   ** makros **
   ************/
  # define __OCL_ASYNC_FUNCTION_OBJECT_HPP__
  
  
  
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
      
      /**********************
       ** member functions **
       **********************/
      virtual
      void
      run_async           () = 0;
      
      virtual
      void
      activate            () = 0;
      
      
    public:
    
      virtual
      bool
      isActive            () = 0;
    
  
  }; // class oclAsyncFunctionObject
  
  
  
# endif // __OCL_ASYNC_FUNCTION_OBJECT_HPP__
