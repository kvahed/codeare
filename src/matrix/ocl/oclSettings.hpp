# ifndef __OCL_SETTINGS_HPP__




  /************
   ** makros **
   ************/

  /* include guard */
  # define __OCL_SETTINGS_HPP__




  /**************
   ** includes **
   **************/
  # include <sstream>




  /**************************
   ** forward declarations **
   **************************/
  
  /* classes */
  class oclError;
  
  /* functions */
  extern
  ostream &
  operator<<        (        ostream &  os,
                      const oclError & err );



  /**********************
   ** type definitions **
   **********************/
  typedef int oclObjectID;

  
  
  
  /****************
   ** enum types **
   ****************/

  enum vclAlgoType
  {
      
    vclSUBTRACT,
    vclMATPROD
  
  };

  
  enum VerbosityLevel
  {
  
    VERB_NONE = 0,
    VERB_LOW,
    VERB_MIDDLE,
    VERB_HIGH
  
  };
  
  enum oclClassID
  {
  
  /* 1. */     OCL_OBSERVED_DATA_OBJECT = 0,
  /* 2. */     OCL_DATA_OBJECT,
  /* 3. */     OCL_DATA_WRAPPER,
  /* 4. */     OCL_CPU_DATA_OBJECT,
  /* 5. */     OCL_GPU_DATA_OBJECT,
  /* 6. */     OCL_TRAITS,
  /* 7. */     OCL_CONNECTION,
  /* 8. */     OCL_FUNCTION_OBSERVER,
  /* 9. */     OCL_FUNCTION_OBJECT,
  /* 10. */    OCL_VIENNACL_OBJECT,
  /* 11. */    OCL_KERNEL_OBJECT,
  /* 12. */    OCL_ASYNC_FUNCTION_OBJECT,
  /* 13. */    OCL_ASYNC_VCL_OBJECT,
  /* 14. */    OCL_ASYNC_KERNEL_OBJECT,
    
    MAX_ID
  
  };




  /*****************
   ** global vars **
   *****************/

  VerbosityLevel verbosity = VERB_NONE; //VERB_HIGH;
  
  static const VerbosityLevel global_verbosity [MAX_ID] =
  {
  
  /* 1. */    VERB_HIGH,
  /* 2. */    VERB_HIGH,
  /* 3. */    VERB_HIGH,
  /* 4. */    VERB_HIGH,
  /* 5. */    VERB_HIGH,
  /* 6. */    VERB_LOW,
  /* 7. */    VERB_HIGH,
  /* 8. */    VERB_HIGH,
  /* 9. */    VERB_HIGH,
  /* 10. */   VERB_HIGH,
  /* 11. */   VERB_HIGH,
  /* 12. */   VERB_HIGH,
  /* 13. */   VERB_HIGH,
  /* 14. */   VERB_HIGH
  
  };

  /* buffer for writing verbosity messages */
  char buffer [100];




  /**************************
   ** function definitions **
   **************************/


  /**
   * @brief               print given message, if v_level matches
   */
  inline
  void
  print_optional          ( const           char *     msg,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << msg << std::endl;
    
    }
  
  }
  
  
  /**
   * @brief               print oclError, if v_level matches
   */
  inline
  void
  print_optional          ( const       oclError &     err,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << err << std::endl;
    
    }
  
  }
  
  
  /**
   * @brief               print given messages, if v_level matches
   */
  inline
  void
  print_optional          ( const           char *    msg1,
                            const           char *    msg2,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << msg1 << msg2 << std::endl;
    
    }
  
  }
  
  
  /**
   * @brief               print given messages, if v_level matches
   */
  inline
  void
  print_optional          ( const           char *    msg1,
                            const           char *    msg2,
                            const           char *    msg3,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << msg1 << msg2 << msg3 << std::endl;
    
    }
  
  }


  /**
   * @brief               print given messages, if v_level matches
   */
  inline
  void
  print_optional          ( const           char *    msg1,
                            const           char *    msg2,
                            const           char *    msg3,
                            const           char *    msg4,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << msg1 << msg2 << msg3 << msg4 << std::endl;
    
    }
  
  }


  /**
   * @brief               print given messages, if v_level matches
   */
  inline
  void
  print_optional          ( const           char *    msg1,
                            const           char *    msg2,
                            const           char *    msg3,
                            const           char *    msg4,
                            const           char *    msg5,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      std::cout << msg1 << msg2 << msg3 << msg4 << msg5 << std::endl;
    
    }
  
  }
  
  
  /**
   * @brief               print given message and one integer
   *
   * @param               msg - !formatted! string
   */
  inline
  void
  print_optional          ( const           char *     msg,
                                             int         i,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      /* create string containing interger */
      sprintf (buffer, msg, i);
    
      std::cout << buffer << std::endl;
    
    }
  
  }
  
  
  
  /**
   * @brief               print given message and two integers
   *
   * @param               msg - !formatted! string
   */
  inline
  void
  print_optional          ( const           char *     msg,
                                             int         i,
                                             int         j,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      /* create string containing interger */
      sprintf (buffer, msg, i, j);
    
      std::cout << buffer << std::endl;
    
    }
  
  }
  


  
# endif // __OCL_SETTINGS_HPP__
