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
      oclKernelObject       (std::string                    kernel_name,
                             oclDataObject  * const * const pp_args,
                             int                            num_args )
                          : oclFunctionObject (pp_args, num_args),
                            m_kernel_name     (kernel_name)
                            
      {
      
        std::cout << "Ctor: \"oclKernelObject\"" << std::endl;
        
        /* TODO */
      }
    
    
      /**
       * @brief             virtual destructor
       */
      virtual
      ~oclKernelObject      ()
      {
      
        std::cout << "Dtor: \"oclKernelObject\"" << std::endl;

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
      
  
  
  }; // class oclKernelObject
  
  
  
  
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
  
    std::cout << "oclKernelObject :: run!" << std::endl;
    
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
    cl::NDRange global_dims (8);
    cl::NDRange local_dims = cl::NullRange;
    oclCon -> runKernel (global_dims, local_dims);
    
    // get data
    for (int i = 0; i < m_num_args; i++)
    {
      mpp_args [i] -> finish ();
    }
    
  }
  
  
  
  
# endif // __OCL_KERNEL_OBJECT_HPP__
