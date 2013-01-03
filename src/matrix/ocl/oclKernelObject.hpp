# ifndef __OCL_KERNEL_OBJECT_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_KERNEL_OBJECT_HPP__




  /**************
   ** includes **
   **************/
  
  // C++ std lib
  # include <string>

  // ocl
  # include "oclSettings.hpp"
  # include "oclDataObject.hpp"
  # include "oclFunctionObject.hpp"

    
  
    
  /****************************
   ** class: oclKernelObject **
   **   (derived)            **
   ****************************/
  class oclKernelObject : public oclFunctionObject
  {
  
    
    
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
   
      std::string           m_kernel_name;
    
    
   
    public:
   
    
      /**
       * @name              Constructors and destructors.
       */
      //@{
      
      
      /**
       * @brief             default constructor
       */
      oclKernelObject       ( const   std::string &               kernel_name,
                                    oclDataObject * const * const     pp_args,
                                              int                    num_args )
                           : oclFunctionObject (pp_args, num_args),
                             m_kernel_name     (kernel_name)
                            
      {
      
        print_optional ("Ctor: \"oclKernelObject\"", v_level);
        
        /* TODO */
        
      }
    
    
      /**
       * @brief             virtual destructor
       */
      virtual
      ~oclKernelObject      ()
      {
      
        print_optional ("Dtor: \"oclKernelObject\"", v_level);

        /* TODO */
        
      }

    
      //@}


      /**
       * @brief             execute kernel
       *
       * @see               defined in oclFunctionObject
       */
      virtual
      void
      run                   ();
      
      

    private:
    
      /* private member for verbosity level of class */
      static const VerbosityLevel v_level;
      
  
  
  }; // class oclKernelObject
  
  
  
  /*************************************
   ** initialize static class members **
   *************************************/
  const VerbosityLevel oclKernelObject :: v_level = global_verbosity [OCL_KERNEL_OBJECT];
  
  
  
  /**************************
   ** function definitions **
   **************************/


  
  /**
   * @brief                 prepare arguments, run kernel, finish arguments
   */
  void
  oclKernelObject ::
  run                       ()
  {
  
    // oclConnection for reuse in this function
    oclConnection * oclCon = oclConnection :: Instance ();
  
    print_optional ("oclKernelObject :: run ()", v_level);
    
    // activate kernel
    oclCon -> activateKernel (m_kernel_name);
    
    // prepare kernel arguments (load to gpu)
    for (int i = 0; i < m_num_args; i++)
    {

      // prepare argument
      mpp_args [i] -> prepare ();
      
      // register argument at kernel
      oclCon -> setKernelArg (i, mpp_args [i]);

    }
    
    // run kernel
    cl::NDRange global_dims (1024);
    cl::NDRange local_dims = cl::NullRange;
    oclCon -> runKernel (global_dims, local_dims);

    // perhaps get data
    for (int i = 0; i < m_num_args; i++)
    {
      mpp_args [i] -> finish ();
    }
    
  }
  
  
  
  
# endif /* __OCL_KERNEL_OBJECT_HPP__ */
