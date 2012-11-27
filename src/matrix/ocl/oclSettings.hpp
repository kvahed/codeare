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
      
    vclSUBTRACT
  
  };

  
  enum VerbosityLevel
  {
  
    VERB_NONE = 0,
    VERB_LOW,
    VERB_MIDDLE,
    VERB_HIGH
  
  };




  /*****************
   ** global vars **
   *****************/

  const VerbosityLevel verbosity = VERB_NONE; //VERB_HIGH;

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
    
      std::cout << msg1 << msg2<< std::endl;
    
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
