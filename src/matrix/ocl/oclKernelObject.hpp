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
  
  
  
  /**************************
   ** forward declarations **
   **************************/
//  class oclDataObject;
//  void oclDataObject :: run ();
  
  
    
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
       * @name      Constructors and destructors.
       */
      //@{
      
      oclKernelObject       (std::string                    kernel_name,
                             oclDataObject  * const * const pp_args,
                             int                            num_args )
                          : oclFunctionObject (pp_args, num_args),
                            m_kernel_name     (kernel_name)
                            
      {
        std::cout << "Ctor: \"oclKernelObject\"" << std::endl;
        /* TODO */
      }
    
      virtual
      ~oclKernelObject ()
      {
        std::cout << "Dtor: \"oclKernelObject\"" << std::endl;
        /* TODO */
      }
    
      //@}
    
      virtual
      void
      run                   ();
      
  
  }; // class oclKernelObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  void
  oclKernelObject ::
  run                       ()
  {
  
    oclConnection * oclCon = oclConnection :: Instance ();
  
    std::cout << "oclKernelObject :: run!" << std::endl;
    
    // activate kernel
    oclCon -> activateKernel (m_kernel_name);
    
    // prepare kernel arguments (load to gpu)
    for (int i = 0; i < m_num_args; i++)
    {
      mpp_args [i] -> prepare (i);
    }
    
    std::cout << "calculate! (Kernel)" << std::endl;
    
    // run kernel
    cl::NDRange global_dims (512);
    cl::NDRange local_dims = cl::NullRange;
    oclCon -> runKernel (global_dims, local_dims);
    
    // get data
    for (int i = 0; i < m_num_args; i++)
    {
      mpp_args [i] -> finish (i);
    }
    
  }
  
  
  
# endif // __OCL_KERNEL_OBJECT_HPP__
