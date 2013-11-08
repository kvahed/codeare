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
   
      std::vector <std::string> m_kernel_names;
      std::vector <cl::Event>   m_events;
      double * mp_time_mem_up;
      double * mp_time_mem_down;
    
    
   
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
                             m_kernel_names    (1, kernel_name),
                             m_events          (),
                             mp_time_mem_up    (new double [1]),
                             mp_time_mem_down  (new double [1]) 
                             
                            
      {
      
        print_optional ("Ctor: \"oclKernelObject\"", v_level);
        
        /* TODO */
        
      }
    
                           
                           
      /**
       * @brief             default constructor
       */
      oclKernelObject       ( const std::vector   <std::string> &               kernel_names,
                                    oclDataObject               * const * const      pp_args,
                                              int                                   num_args )
                           : oclFunctionObject (pp_args, num_args),
                             m_kernel_names    (kernel_names),
                             m_events          (),
                             mp_time_mem_up    (new double [kernel_names.size ()]),
                             mp_time_mem_down  (new double [kernel_names.size ()])
                            
      {
      
        print_optional ("Ctor: \"oclKernelObject\"", v_level);
        
        /* TODO */
        
      }
      
      
      /**
       * @brief             copy constructor
       */
      oclKernelObject       (const oclKernelObject & cp_obj)
                            : oclFunctionObject (cp_obj.mpp_args, m_num_args),
                              m_kernel_names    (cp_obj.m_kernel_names),
                              m_events          (),
                              mp_time_mem_up    (new double [cp_obj.m_kernel_names.size ()]),
                              mp_time_mem_down  (new double [cp_obj.m_kernel_names.size ()])
      {
        
        print_optional ("Ctor: \"oclKernelObject\"", v_level);
        
        // copy measured memory transfer times
        for (int i = 0; i < m_kernel_names.size (); ++i)
        {
          mp_time_mem_up [i] = cp_obj.mp_time_mem_up [i];
          mp_time_mem_down [i] = cp_obj.mp_time_mem_down [i];
        }
        
      }
      
      
                           
      /**
       * @brief             virtual destructor
       */
      virtual
      ~oclKernelObject      ()
      {
      
        print_optional ("Dtor: \"oclKernelObject\"", v_level);

        /* TODO */
        delete [] mp_time_mem_up;
        delete [] mp_time_mem_down;
        
      }

    
      //@}
      
      
      /**
       * @name              Operators.
       */
      //@{
      
      oclKernelObject &
      operator=             (const oclKernelObject & assign_obj)
      {
        this -> m_events = assign_obj.m_events;
        this -> m_kernel_names = assign_obj.m_kernel_names;
        this -> m_num_args = assign_obj.m_num_args;
        this -> mpp_args = assign_obj.mpp_args;
        delete this -> mp_time_mem_up;
        delete this -> mp_time_mem_down;
        this -> mp_time_mem_up = new double [assign_obj.m_kernel_names.size ()];
        this -> mp_time_mem_down = new double [assign_obj.m_kernel_names.size ()];
        for (int i = 0; i < assign_obj.m_kernel_names.size (); ++i)
        {
          this -> mp_time_mem_up [i] = assign_obj.mp_time_mem_up [i];
          this -> mp_time_mem_down [i] = assign_obj.mp_time_mem_down [i];
        }
        return *this;
      }
      
      //@}


      /**
       * @brief             execute kernel
       *
       * @see               defined in oclFunctionObject
       */
      virtual
      void
      run                   (const LaunchInformation &);
      
      
      /**
       * @brief             Retrieve profiling information on kernel execution.
       */
      virtual
      const ProfilingInformation
      getProfilingInformation (const int i)
      const;
      
      

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
  run                       (const LaunchInformation & lc)
  {
  
    // oclConnection for reuse in this function
    oclConnection * oclCon = oclConnection :: Instance ();
    
    double * p_mem_time = mp_time_mem_up;
    
    for (std::vector<std::string>::const_iterator it_kernel_name = m_kernel_names.begin (); it_kernel_name != m_kernel_names.end (); it_kernel_name++, p_mem_time++)
    {

      print_optional("oclKernelObject :: run ( \"", it_kernel_name->c_str(), "\" )", v_level);

      // activate kernel
      oclCon -> activateKernel(*it_kernel_name);

      // prepare kernel arguments (load to gpu)
      for (int i = 0; i < m_num_args; i++)
      {
        
        // prepare argument
        *p_mem_time += mpp_args [i] -> prepare();
        
        // register argument at kernel
        oclCon -> setKernelArg(i, mpp_args [i]);

      }

      // run kernel
      cl::NDRange global_dims(lc.global_x, lc.global_y);
      cl::NDRange local_dims(lc.local_x, lc.local_y);
      cl::Event event = oclCon -> runKernel(global_dims, local_dims);

      m_events.push_back (event);
      
    }
    
    for (int i = 0; i < m_events.size (); i++)
    {
      try {
      m_events [i].wait ();
      } catch (cl::Error & cle)
      {
        std::cerr << oclError (cle.err (), " oclKernelObject::run ") << std::endl;
      }
    }
    
    // perhaps get data
    for (int i = 0; i < m_num_args; i++)
    {
      double time_mem = mpp_args [i] -> finish ();
      for (int j = 0; j < m_events.size (); ++j)
        mp_time_mem_down [j] += time_mem;
    }
        
  }
  
  
  const ProfilingInformation
  oclKernelObject ::
  getProfilingInformation     (const int i)
  const
  {
    ProfilingInformation pi = oclConnection :: Instance () -> getProfilingInformation (m_events [i]);
    pi.time_mem_up = mp_time_mem_up [i];
    pi.time_mem_down = mp_time_mem_down [i];
    return pi;
  }
  
  
  
  
# endif /* __OCL_KERNEL_OBJECT_HPP__ */
