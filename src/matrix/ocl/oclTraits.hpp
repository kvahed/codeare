# ifndef __OCL_TRAITS_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_TRAITS_HPP__
  
  

  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclGPUDataObject.hpp"
    
  

  
  /***********************
   ** struct: oclTraits **
   **   (base struct)   **
   ***********************/
  template <class T>
  struct oclTraits
  {

    /* -- */

  }; // struct oclTraits

  
  
  
  /*******************************
   ** struct: oclTraits <float> **
   **   (single precision)      **
   *******************************/
  template <>
  struct oclTraits <float>
  {
    

    
    /**********************
     ** type definitions **
     **********************/
    typedef float elem_type;

    
    
    /**
     * @brief                       Create oclGPUDataObject.
     */
    static inline
    oclDataWrapper <elem_type> *
    make_GPU_Obj                    (elem_type * const cpu_arg, const size_t & num_elems)
    {
    
      std::cout << "make_GPU_Obj <float> (create new)" << std::endl;
    
      return new oclGPUDataObject <elem_type> (cpu_arg, num_elems * sizeof (elem_type));
      
    }
    
    
    /**
     * @brief                      Create oclGPUDataObject
     *                             with same state as given oclDataWrapper.
     */
    static inline
    oclDataWrapper <elem_type> *
    make_GPU_Obj                   (elem_type * const cpu_arg, const size_t & num_elems, const oclDataWrapper <elem_type> & obj)
    {
    
      std::cout << "make_GPU_Obj <float> (copy obj's state)" << std::endl;
      
      /**
       * local vars
       */
      size_t size = num_elems * sizeof (elem_type);
      
      // check if sizes fit (for !some! more control (or safety))
      if (size != obj.getSize ())
      {
      
        return NULL;
      
      }
      else
      {
      
        return new oclGPUDataObject <elem_type> (cpu_arg, size, obj);
      
      }
    
    }
    
    
    
    /**
     * @name                        operators
     */
    //@{
    
    
    /**
     * @brief                       elementwise addition
     */
    static inline
    const oclError &
    ocl_operator_add                (oclDataObject * const arg1, oclDataObject * const arg2, oclDataObject * const sum, size_t size)
    {
      
      std::cout << "oclTraits <float> :: ocl_operator_add" << std::endl;
            
      // create array of function arguments
      oclDataObject ** args = (oclDataObject **) malloc (4 * sizeof (oclDataObject *));
      args [0] = arg1;
      args [1] = arg2;
      args [2] = sum;
      args [3] = new oclGPUDataObject <size_t> (& size, 1 * sizeof (size_t));

      // create function object
      oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject ("add", args, 4, oclConnection::KERNEL, oclConnection::SYNC);
      
      // execute function object
      op_obj -> run ();

      // clear memory      
      delete op_obj;      
      delete args [3];
      free (args);
                  
    }
    
    
    /**
     * @brief                       elementwise subtraction
     */
    static inline
    const oclError &
    ocl_operator_subtract           (oclDataObject * const arg1, oclDataObject * const arg2, oclDataObject * const diff, size_t size)
    {
      
      std::cout << "oclTraits <float> :: ocl_operator_subtract" << std::endl;
      
      // create array of function arguments
      oclDataObject ** args = (oclDataObject **) malloc (4 * sizeof (oclDataObject *));
      args [0] = arg1;
      args [1] = arg2;
      args [2] = diff;
      args [3] = new oclGPUDataObject <size_t> (& size, 1 * sizeof (size_t));
            
      // create function object
      oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject ("subtract", args, 4, oclConnection::VCL, oclConnection::SYNC);
      
      // execute function object
      op_obj -> run ();
      
      // clear memory
      delete op_obj;
      delete args [3];
      free (args);
      
    }
    
    
    //@}
    
    
    
  };
  
  
  
  
  /********************************
   ** struct: oclTraits <size_t> **
   **   (size_t)                 **
   ********************************/
  template <>
  struct oclTraits <size_t>
  {
    
    
    
    /**********************
     ** type definitions **
     **********************/
    typedef size_t elem_type;
    
    
    
    /**
     * @brief                       Create oclDataObject
     */
    static inline
    oclDataWrapper <elem_type> *
    make_GPU_Obj                    (elem_type * const cpu_arg, const size_t & num_elems)
    {
    
      std::cout << "make_GPU_Obj <size_t>" << std::endl;
    
      return new oclGPUDataObject <elem_type> (cpu_arg, num_elems * sizeof (elem_type));
      
    }
    
    
    
  };
  
  
  
  
# endif // __OCL_TRAITS_HPP__
