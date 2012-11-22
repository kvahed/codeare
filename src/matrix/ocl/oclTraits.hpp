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
  # include "oclSettings.hpp"
      
  

  
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
    

    private:

    
      /**********************
       ** type definitions **
       **********************/
      typedef float elem_type;


      /**
       * @name                        basic operator algos
       */
      //@{
    
    
      /**
       * @brief                       execute specified kernel with 3 arguments
       */
      static inline
      const oclError &
      ocl_basic_operator_kernel_3              ( const          char * const kernel_name,
                                                       oclDataObject * const        arg1,
                                                       oclDataObject * const        arg2,
                                                       oclDataObject * const      result,
                                                                 int           num_elems )
      {
    
        // number of kernel arguments
        const int num_args = 4;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = arg2;
        args [2] = result;
        args [3] = new oclGPUDataObject <int> (& num_elems, 1);

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (kernel_name, args, num_args, oclConnection::KERNEL, oclConnection::SYNC);
      
        // execute function object
        op_obj -> run ();

        // clear memory      
        delete op_obj;      
        delete args [3];
        free (args);
    
      }
    
    
      /**
       * @brief                      execute specified kernel with 2 arguments
       */
      static inline
      const oclError &
      ocl_basic_operator_kernel_2              ( const          char * const kernel_name,
                                                       oclDataObject * const        arg1,
                                                       oclDataObject * const        arg2,
                                                                 int           num_elems )
      {
    
        // number of kernel arguments
        const int num_args = 3;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = arg2;
        args [2] = new oclGPUDataObject <int> (& num_elems, 1);

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (kernel_name, args, num_args, oclConnection::KERNEL, oclConnection::SYNC);
      
        // execute function object
        op_obj -> run ();

        // clear memory      
        delete op_obj;      
        delete args [3];
        free (args);
    
      }
    
    
      /**
       * @brief                       execute specified ViennaCl algorithm with 3 arguments
       */
      static inline
      const oclError &
      ocl_basic_operator_vclAlgo_3             ( const   vclAlgoType            vcl_algo,
                                                       oclDataObject * const        arg1,
                                                       oclDataObject * const        arg2,
                                                       oclDataObject * const      result,
                                                                 int           num_elems )
      {
    
        // number of kernel arguments
        const int num_args = 4;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = arg2;
        args [2] = result;
        args [3] = new oclGPUDataObject <int> (& num_elems, 1);

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (vcl_algo, args, num_args, oclConnection::VCL, oclConnection::SYNC);
      
        // execute function object
        op_obj -> run ();

        // clear memory      
        delete op_obj;      
        delete args [3];
        free (args);
    
      }        
    
    
      //@}


    public:
    
    
      /**
       * @name                        memory management
       */
      //@{
    
    
      /**
       * @brief                       Create oclGPUDataObject.
       */
      static inline
      oclDataWrapper <elem_type> *
      make_GPU_Obj                    (      elem_type * const   cpu_arg,
                                       const    size_t &       num_elems)
      {
    
        std::cout << "make_GPU_Obj <float> (create new)" << std::endl;
    
        return new oclGPUDataObject <elem_type> (cpu_arg, num_elems);
      
      }
    
    
      /**
       * @brief                      Create oclGPUDataObject
       *                             with same state as given oclDataWrapper.
       */
      static inline
      oclDataWrapper <elem_type> *
      make_GPU_Obj                   (           elem_type             *  const     cpu_arg,
                                      const         size_t             &          num_elems,
                                            oclDataWrapper <elem_type> &                obj,
                                                      bool                      keep_buffer = false)
      {
      
        std::cout << "make_GPU_Obj <float> (copy obj's state)" << std::endl;
      
        // check if sizes fit (for !some! more control (or safety))
        if (num_elems != obj.getNumElems ())
        {
      
          return NULL;
      
        }
        else
        {
      
          // create oclDataObject to return
          oclDataWrapper <elem_type> * cp_obj = new oclGPUDataObject <elem_type> (cpu_arg, num_elems, obj, keep_buffer);
        
          /* update GPU buffer with data if possible */
          if (keep_buffer && obj.bufferCopyable ())
          {
       
            // number of kernel arguments
            const int num_args = 3;
    
            int tmp_num_elems = num_elems;
    
            // create array of function arguments
            oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
            args [0] = cp_obj;
            args [1] = & obj;
            args [2] = new oclGPUDataObject <int> (& tmp_num_elems, 1);

            // create function object
            oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> ("copy_buffer", args, num_args, oclConnection::KERNEL, oclConnection::SYNC);
      
            // execute function object
            op_obj -> run ();

            // clear memory      
            delete op_obj;      
            delete args [2];
            free (args);
  
          }
        
          return cp_obj;
      
        }
    
      }
    
    
      //@}
    
    
      /**
       * @name                        operators
       */
      //@{
    
    
      /**
       * @brief                       elementwise addition
       */
      static inline
      const oclError &
      ocl_operator_add                ( oclDataObject * const      arg1,
                                        oclDataObject * const      arg2,
                                        oclDataObject * const       sum,
                                                  int         num_elems )
      {
      
        std::cout << "oclTraits <float> :: ocl_operator_add" << std::endl;
      
        ocl_basic_operator_kernel_3 ("add", arg1, arg2, sum, num_elems);
      
      }
        
    
      /**
       * @brief                       elementwise subtraction
       */
      static inline
      const oclError &
      ocl_operator_subtract           (oclDataObject * const      arg1,
                                       oclDataObject * const      arg2,
                                       oclDataObject * const      diff,
                                                 int         num_elems)
      {
      
        std::cout << "oclTraits <float> :: ocl_operator_subtract" << std::endl;
      
        ocl_basic_operator_vclAlgo_3 (vclSUBTRACT, arg1, arg2, diff, num_elems);
      
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
