# ifndef __OCL_SETTINGS_HPP__




  /************
   ** makros **
   ************/

  /* include guard */
  # define __OCL_SETTINGS_HPP__


//  # define __USE_VIENNA_CL__




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
# ifdef __USE_VIENNA_CL__
  enum vclAlgoType
  {

    vclSUBTRACT,
    vclMATPROD
  
  };
# endif
  
  enum oclAMDBlasType
  {
  
    amdblasGEMM,
    amdblasGEMV
  
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
  
  
  // base path for kernel files
  static std::string base_kernel_path = "/localdata/djoergens/projects/CoDEARE/";
  
  
  
  /***********************
   ** class definitions **
   ***********************/

  
  /**
   * @brief  Struct containing profiling information.
   */
  struct ProfilingInformation
  {
    double time_start;
    double time_end;
    double time_mem_up;
    double time_mem_down;
  };
  
  
  
  /**
   * @brief  Struct containing information on kernel launches.
   */
  struct LaunchInformation
  {
    int local_x;
    int local_y;
    int global_x;
    int global_y;
    LaunchInformation (const int loc_x, const int loc_y, const int glob_x, const int glob_y)
    {
      local_x = loc_x;
      local_y = loc_y;
      global_x = glob_x;
      global_y = glob_y;
    }
    LaunchInformation (const LaunchInformation & lc)
    {
      local_x = lc.local_x;
      local_y = lc.local_y;
      global_x = lc.global_x;
      global_y = lc.global_y;
    }
    bool
    operator==        (const LaunchInformation & lc)
    const
    {
      return this -> local_x == lc.local_x
              && this -> local_y == lc.local_y
              && this -> global_x == lc.global_x
              && this -> global_y == lc.global_y;
    }
  };
  
  
  struct PerformanceInformation
  {
    std::string kernel_name;
    LaunchInformation lc;
    std::string information;
    double time_exec;
    double time_mem_up;
    double time_mem_down;
    double parameter;
    PerformanceInformation (const std::string & k_name, const LaunchInformation & _lc, const std::string & inf, const double & t_exec, const double & t_mem_up, const double & t_mem_down, const double & param)
      : lc (_lc)
    {
      kernel_name = k_name;
      information = inf;
      time_exec = t_exec;
      time_mem_up = t_mem_up;
      time_mem_down = t_mem_down;
      parameter = param;
    }
    PerformanceInformation &
    operator+=      (const PerformanceInformation & pi)
    {
      if (0 == this -> kernel_name.compare (pi.kernel_name) && this -> lc == pi.lc && this -> information.compare (pi.information) == 0)
      {
        this -> time_exec = (this -> time_exec+pi.time_exec)/2;
        this -> time_mem_up = (this -> time_mem_up+pi.time_mem_up)/2;
        this -> time_mem_down = (this -> time_mem_down+pi.time_mem_down)/2;
        this -> parameter = (this -> parameter+pi.parameter)/2;
      }
      else
        throw -1;
      return *this;
    }
  };
  
  
  PerformanceInformation
  operator+         (const PerformanceInformation & pi1, const PerformanceInformation & pi2)
  {
    if (0 == pi1.kernel_name.compare (pi2.kernel_name) && pi1.lc == pi2.lc && pi1.information.compare (pi2.information) == 0)
      return PerformanceInformation (pi1.kernel_name, pi1.lc, pi1.information, (pi1.time_exec+pi2.time_exec)/2, (pi1.time_mem_up+pi2.time_mem_up)/2, (pi1.time_mem_down+pi2.time_mem_down)/2, (pi1.parameter+pi2.parameter)/2);
    else
      throw -1;
  }
  
  
  std::ostream &
  operator<<        (        std::ostream &  os,
                      const PerformanceInformation & pi )
  {
    os << " **> Kernel: " << pi.kernel_name << " <**" << std::endl;
    os << "   local size: " << pi.lc.local_x << " x " << pi.lc.local_y << std::endl;
    os << "   global size: " << pi.lc.global_x << " x " << pi.lc.global_y << std::endl;
    os << "   Execution time in seconds: " << pi.time_exec << " s " << std::endl;
    os << "   Memory transfer time in seconds: " << (pi.time_mem_up + pi.time_mem_down) << " s " << std::endl;
    os << "  " << pi.information << ": " << pi.parameter << std::endl;
    return os;
  }
  
  

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
  
  
  
    /**
   * @brief               print given message and three integers
   *
   * @param               msg - !formatted! string
   */
  inline
  void
  print_optional          ( const           char *     msg,
                                             int         i,
                                             int         j,
                                             int         k,
                                  VerbosityLevel   v_level )
  {
  
    /* check verbosity level */
    if (v_level <= verbosity)
    {
    
      /* create string containing interger */
      sprintf (buffer, msg, i, j, k);
    
      std::cout << buffer << std::endl;
    
    }
  
  }
  


  
# endif // __OCL_SETTINGS_HPP__
