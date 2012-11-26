# ifndef __OCL_ERROR_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_ERROR_HPP__
  
  
  
  
  /**************
   ** includes **
   **************/
  # include <sstream>
  
  
  
  
  /**************************
   ** function definitions **
   **************************/

  /**
   * @brief               translate OpenCl errors
   */
  const char*
  getErrorString             (cl_int e)
  {

    static const char* errorString[] = {
      "CL_SUCCESS",
      "CL_DEVICE_NOT_FOUND",
      "CL_DEVICE_NOT_AVAILABLE",
      "CL_COMPILER_NOT_AVAILABLE",
      "CL_MEM_OBJECT_ALLOCATION_FAILURE",
      "CL_OUT_OF_RESOURCES",
      "CL_OUT_OF_HOST_MEMORY",
      "CL_PROFILING_INFO_NOT_AVAILABLE",
      "CL_MEM_COPY_OVERLAP",
      "CL_IMAGE_FORMAT_MISMATCH",
      "CL_IMAGE_FORMAT_NOT_SUPPORTED",
      "CL_BUILD_PROGRAM_FAILURE",
      "CL_MAP_FAILURE",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "CL_INVALID_VALUE",
      "CL_INVALID_DEVICE_TYPE",
      "CL_INVALID_PLATFORM",
      "CL_INVALID_DEVICE",
      "CL_INVALID_CONTEXT",
      "CL_INVALID_QUEUE_PROPERTIES",
      "CL_INVALID_COMMAND_QUEUE",
      "CL_INVALID_HOST_PTR",
      "CL_INVALID_MEM_OBJECT",
      "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
      "CL_INVALID_IMAGE_SIZE",
      "CL_INVALID_SAMPLER",
      "CL_INVALID_BINARY",
      "CL_INVALID_BUILD_OPTIONS",
      "CL_INVALID_PROGRAM",
      "CL_INVALID_PROGRAM_EXECUTABLE",
      "CL_INVALID_KERNEL_NAME",
      "CL_INVALID_KERNEL_DEFINITION",
      "CL_INVALID_KERNEL",
      "CL_INVALID_ARG_INDEX",
      "CL_INVALID_ARG_VALUE",
      "CL_INVALID_ARG_SIZE",
      "CL_INVALID_KERNEL_ARGS",
      "CL_INVALID_WORK_DIMENSION",
      "CL_INVALID_WORK_GROUP_SIZE",
      "CL_INVALID_WORK_ITEM_SIZE",
      "CL_INVALID_GLOBAL_OFFSET",
      "CL_INVALID_EVENT_WAIT_LIST",
      "CL_INVALID_EVENT",
      "CL_INVALID_OPERATION",
      "CL_INVALID_GL_OBJECT",
      "CL_INVALID_BUFFER_SIZE",
      "CL_INVALID_MIP_LEVEL",
      "CL_INVALID_GLOBAL_WORK_SIZE",
    };

    const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

    const int index = -e;

    return (index >= 0 && index < errorCount) ? errorString[index] : "";

  }
  
  
  
  
  /**********************
   ** class definition **
   **********************/
  
  /**
   * @brief           wrapper for errors, handling errors
   */
  class oclError
  {


    public:
    
      /**
       * @brief       enum for error types
       */
      enum ErrorType
      {
      
        OPENCL_ERROR,
        CUSTOM_ERROR,
        WARNING
      
      };
    
    
      /**
       * @name        constructors and destructors
       */
      //@{
    
    
      /**
       * @brief       constructor, own error message
       */
      oclError        ( const      char * const    msg,
                        const      char * const source,
                        const ErrorType           type = CUSTOM_ERROR )
                     : m_msg                                                     (msg),
                       m_source (((std::string ()).append ("\n -> ")).append (source)),
                       m_type                                                   (type)
      {
      
        print_optional ("CTOR: \"oclError\"", VERB_HIGH);
        print_optional (" *!* Error: ", msg, VERB_HIGH);
      
      }
      
      
      /**
       * @brief       constructor, own error message
       */
      oclError        ( const std::string      msg,
                        const        char * source,
                        const   ErrorType     type = CUSTOM_ERROR )
                     : m_msg                                                     (msg),
                       m_source (((std::string ()).append ("\n -> ")).append (source)),
                       m_type                                                   (type)
      {
      
        print_optional ("CTOR: \"oclError\"", VERB_HIGH);
        print_optional (" *!* Error: ", msg.c_str (), VERB_HIGH);
      
      }
      
      
      /**
       * @brief       constructor, OpenCl error code
       */
      oclError        ( const cl_int     code,
                        const   char * source )
                     : m_msg                                                  (getErrorString (code)),
                       m_source                (((std::string ()).append ("\n -> ")).append (source)),
                       m_type                                                          (OPENCL_ERROR)
      {
      
        print_optional ("CTOR: \"oclError\"", VERB_HIGH);
        print_optional (" *!* OpenCl error: ", m_msg.c_str (), VERB_HIGH);
      
      }
      
      
      /**
       * @brief       constructor, wrap oclError
       */
      oclError        ( const oclError & err,
                        const     char * source )
                     : m_msg                                                           (err.m_msg),
                       m_source (((std::string (err.m_source)).append ("\n -> ")).append (source)),
                       m_type                                                         (err.m_type)
      {
      
        print_optional ("CTOR: \"oclError\"", VERB_HIGH);
        print_optional (" *!* Error: ", m_msg.c_str (), VERB_HIGH);
      
      }
      
      
      //@}
    
    
      /*************************
       ** friend declarations **
       *************************/
      friend
      ostream &
      operator<<      (        ostream &  os,
                        const oclError & err );
    
    
    private:
    
      /**********************
       ** member variables **
       **********************/
      std::string     m_msg;
      std::string     m_source;
      ErrorType       m_type;


  };
  
  
  
  
  /*********************************
   ** global function definitions **
   *********************************/

  /**
   * @brief           output operator for oclError
   */
  ostream &
  operator<<          (        ostream &  os,
                        const oclError & err )
  {
  
    switch (err.m_type)
    {
    
      case oclError :: OPENCL_ERROR:
      
        /* print OpenCl error message */
        os << " *!* OpenCl error: \"";
  
        break;
        
      case oclError :: CUSTOM_ERROR:
      
        /* print custom error message */
        os << " *!* Custom error: \"";
  
        break;

      case oclError :: WARNING:
      
        /* print warning */
        os << " *!* Warning: \"";
        
        break;
        
    }
    
    /* append particular message */
    os << err.m_msg << "\" *!* ";
    
    /* print trace */
    os << err.m_source << std::endl;
    
    return os;
  
  }
  
  
  
  
# endif /* __OCL_ERROR_HPP__ */
