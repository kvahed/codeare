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
      
  

  /******************************
   ** struct: elem_type_traits **
   **     (base struct)        **
   ******************************/
  template <class T>
  struct elem_type_traits
  {
    
    /* -- */
    
  }; // struct elem_type_traits <T>


  /******************************
   ** struct: elem_type_traits **
   **     (spec: float)        **
   ******************************/
  template <>
  struct elem_type_traits <float>
  {
  
    public:

      typedef float elem_type;
      
      static inline
      const char *
      print_elem_type       ( )
      {
        return "float";
      }
  
  }; // struct elem_type_traits <float>


  /******************************
   ** struct: elem_type_traits **
   **     (spec: float)        **
   ******************************/
  template <>
  struct elem_type_traits <double>
  {
  
    public:

      typedef double elem_type;
      
      static inline
      const char *
      print_elem_type       ( )
      {
        return "double";
      }
  
  }; // struct elem_type_traits <double>


  /******************************
   ** struct: elem_type_traits **
   **     (spec: size_t)       **
   ******************************/
  template <>
  struct elem_type_traits <size_t>
  {
  
    public:

      typedef size_t elem_type;
      
      static inline
      const char *
      print_elem_type       ( )
      {
        return "size_t";
      }
  
  }; // struct elem_type_traits <size_t>
  
  
  template <>
  struct elem_type_traits <bool>
  {
  
    public:

      typedef bool elem_type;
      
      static inline
      const char *
      print_elem_type       ( )
      {
        return "bool";
      }
  
  }; // struct elem_type_traits <bool>


  
  /***************************
   ** struct: oclOperations **
   **   (base struct)       **
   ***************************/
  template <class T, class trait = elem_type_traits <T> >
  struct oclOperations
  {
    

    private:

    
      /**********************
       ** type definitions **
       **********************/
      typedef typename trait :: elem_type elem_type;
      
      
      /*********************
       ** local variables **
       *********************/
      
      /* verbosity of operators */
      static
      const VerbosityLevel op_v_level = VERB_LOW;


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
        
        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {
        
          stringstream msg;
          msg << "oclOperations <" << trait :: print_elem_type () << "> :: ocl_basic_operator_kernel_3";
        
          throw oclError (oclError (err, msg.str ().c_str ()), kernel_name);
        
        }

        // clear memory
        delete op_obj;      
        delete args [3];
        free (args);
    
      }
    
    
      /**
       * @brief                      execute specified kernel with 1 arguments and 1 scalar
       */
      static inline
      const oclError &
      ocl_basic_operator_kernel_11              ( const          char * const kernel_name,
                                                        oclDataObject * const        arg1,
                                                            elem_type                arg2,
                                                                  int           num_elems )
      {
    
        // number of kernel arguments
        const int num_args = 3;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = new oclGPUDataObject <elem_type> (&      arg2, 1);
        args [2] = new oclGPUDataObject       <int> (& num_elems, 1);

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (kernel_name, args, num_args, oclConnection::KERNEL, oclConnection::SYNC);
        
        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {
        
          stringstream msg;
          msg << "oclTraits <" << trait :: print_elem_type () << "> :: ocl_basic_operator_kernel_11";
        
          throw oclError (oclError (err, msg.str ().c_str ()), kernel_name);
        
        }

        // clear memory      
        delete op_obj;      
        delete args [1];
        delete args [2];
        free (args);
    
      }
    
    
      /**
       * @brief                      execute specified kernel with 2 arguments and 1 scalar
       */
      static inline
      const oclError &
      ocl_basic_operator_kernel_21              ( const          char * const kernel_name,
                                                        oclDataObject * const        arg1,
                                                            elem_type                arg2,
                                                        oclDataObject * const        arg3,
                                                                  int           num_elems )
      {
    
        // number of kernel arguments
        const int num_args = 4;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = new oclGPUDataObject <elem_type> (&      arg2, 1);
        args [2] = arg3;
        args [3] = new oclGPUDataObject       <int> (& num_elems, 1);

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (kernel_name, args, num_args, oclConnection::KERNEL, oclConnection::SYNC);
        
        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {
        
          stringstream msg;
          msg << "oclTraits <" << trait :: print_elem_type () << "> :: ocl_basic_operator_kernel_21";
        
          throw oclError (oclError (err, msg.str ().c_str ()), kernel_name);
        
        }

        // clear memory      
        delete op_obj;      
        delete args [1];
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
      
        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {
        
          stringstream msg;
          msg << "oclTraits <" << trait :: print_elem_type () << "> :: ocl_basic_operator_kernel_2";
        
          throw oclError (oclError (err, msg.str ().c_str ()), kernel_name);
        
        }

        // clear memory      
        delete op_obj;      
        delete args [2];
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
      
        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {

          stringstream msg;
          msg << "oclTraits <" << trait :: print_elem_type () << "> :: ocl_basic_operator_vclAlgo_3";
        
          throw oclError (err, msg.str ().c_str ());
        
        }

        // clear memory      
        delete op_obj;      
        delete args [3];
        free (args);
    
      }
      
      
      /**
       * @brief                       execute specified ViennaCl algorithm with 3 arguments + 5 scalars
       */
      static inline
      const oclError &
      ocl_basic_operator_vclAlgo_35            ( const   vclAlgoType            vcl_algo,
                                                       oclDataObject * const        arg1,
                                                       oclDataObject * const        arg2,
                                                       oclDataObject * const      result,
                                                                 int                  s1,
                                                                 int                  s2,
                                                                 int                  s3,
                                                                 int                  s4,
                                                                 int                  s5 )
      {
    
        // number of kernel arguments
        const int num_args = 3;
        const int num_scalars = 5;
    
        // create array of function arguments
        oclDataObject ** args = (oclDataObject **) malloc (num_args * sizeof (oclDataObject *));
        args [0] = arg1;
        args [1] = arg2;
        args [2] = result;

        // create array of scalars
        int * scalars = (int *) malloc (3 * sizeof (oclDataObject *));
        scalars [0] = s1;
        scalars [1] = s2;
        scalars [2] = s3;
        scalars [3] = s4;
        scalars [4] = s5;

        // create function object
        oclFunctionObject * op_obj = oclConnection :: Instance () -> makeFunctionObject <elem_type> (vcl_algo, args, num_args, oclConnection::VCL, oclConnection::SYNC, num_scalars, scalars);

        try
        {

          // activate precision mode for type elem_type
          oclConnection :: Instance () -> activate <elem_type> ();

          // execute function object
          op_obj -> run ( );

        }
        catch (oclError & err)
        {
        
          stringstream msg;
          msg << "oclTraits <" << trait :: print_elem_type () << "> :: ocl_basic_operator_vclAlgo_35";
          
          throw oclError (err, msg.str ().c_str ());
        
        }

        // clear memory      
        delete op_obj;
        free (scalars);
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
    
        print_optional ("make_GPU_Obj <", trait :: print_elem_type (), "> (create new)", VERB_HIGH);
    
        return new oclGPUDataObject <elem_type> (cpu_arg, num_elems);
      
      }
    
    
      /**
       * @brief                      Create oclGPUDataObject
       *                             with same state as given oclDataWrapper.
       */
      static inline
      oclDataWrapper <elem_type> *
      make_GPU_Obj                   (                      elem_type             *  const   cpu_arg,
                                      const                    size_t             &        num_elems,
                                                       oclDataWrapper <elem_type> &              obj,
                                            oclDataObject :: CopyMode                      copy_mode = oclDataObject :: NO_BUFFER)
      {
      
        print_optional ("make_GPU_Obj <", trait :: print_elem_type (), "> (copy obj's state)", VERB_HIGH);
      
        // check if sizes fit (for !some! more control (or safety))
        if (num_elems != obj.getNumElems ())
        {
      
          return NULL;
      
        }
        else
        {
      
          // create oclDataObject to return
          oclDataWrapper <elem_type> * cp_obj = new oclGPUDataObject <elem_type> (cpu_arg, num_elems, obj, copy_mode);
        
          /* update GPU buffer with data if possible */
          if (copy_mode == oclDataObject :: COPY_BUFFER  &&  obj.bufferCopyable ())
          {
  
            ocl_basic_operator_kernel_2 ("copy_buffer", cp_obj, & obj, num_elems);
  
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
       * @brief                       Matrix product.
       *
       * @param  arg1                 Address of first factor   (m x k matrix).
       * @param  arg2                 Address of second factor  (k x n matrix).
       * @param  prod                 Address of product        (m x n matrix).
       * @param  m                    First dimension of product.
       * @param  k                    Inner dimension.
       * @param  n                    Second dimension of product.
       * @param  trans1               1 -> Transpose first matrix.
       * @param  trans2               1 -> Transpose second matrix.
       */
      static inline
      const oclError &
      ocl_operator_matprod            ( oclDataObject * const     arg1,
                                        oclDataObject * const     arg2,
                                        oclDataObject * const     prod,
                                                  int                m,
                                                  int                k,
                                                  int                n,
                                                  int           transA,
                                                  int           transB )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_matprod", op_v_level);
        ocl_basic_operator_vclAlgo_35 (vclMATPROD, arg1, arg2, prod, m, k, n, transA, transB);
      
      }
    
    
      /**
       * @brief                       Elementwise addition of two vectors.
       *
       * @param  arg1                 Address of first ocl data object.
       * @param  arg2                 Address of second ocl data object.
       * @param  sum                  Address of resulting ocl data object.
       * @param  num_elems            Number of elements.
       */
      static inline
      const oclError &
      ocl_operator_add                ( oclDataObject * const      arg1,
                                        oclDataObject * const      arg2,
                                        oclDataObject * const       sum,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_add", op_v_level);
        ocl_basic_operator_kernel_3 ("add", arg1, arg2, sum, num_elems);
        
      }
        
    
      /**
       * @brief                       Elementwise subtraction of two vectors.
       *
       * @param  arg1                 Address of first ocl data object.
       * @param  arg2                 Address of second ocl data object.
       * @param  diff                 Address of resulting ocl data object.
       * @param  num_elems            Number of elements.
       */
      static inline
      const oclError &
      ocl_operator_subtract           ( oclDataObject * const      arg1,
                                        oclDataObject * const      arg2,
                                        oclDataObject * const      diff,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_subtract", op_v_level);
        ocl_basic_operator_vclAlgo_3 (vclSUBTRACT, arg1, arg2, diff, num_elems);
      
      }
      
      
      /**
       * @brief                       Elementwise increment of vector.
       *
       * @param  arg1                 Address of vector's data object to be incremented.
       * @param  inc                  Scalar (increment).
       * @param  num_elems            Number of vector's elements.
       */
      static inline
      const oclError &
      ocl_operator_inc                ( oclDataObject * const      arg1,
                                            elem_type               inc,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_inc", op_v_level);
        ocl_basic_operator_kernel_11 ("inc", arg1, inc, num_elems);
      
      }


      /**
       * @brief                       Elementwise decrement of vector.
       *
       * @param  arg1                 Address of vector's data object to be decremented.
       * @param  dec                  Scalar (decrement).
       * @param  num_elems            Number of vector's elements.
       */
      static inline
      const oclError &
      ocl_operator_dec                ( oclDataObject * const      arg1,
                                            elem_type               dec,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_dec", op_v_level);
        
        /* use increment kernel with inverse decrement */
        ocl_basic_operator_kernel_11 ("inc", arg1, -1 * dec, num_elems);
      
      }
      
      
      /**
       * @brief                       Scalar equality.
       *
       * @param  arg1                 Address of vector's data object.
       * @param  scalar               Scalar.
       * @param  result               Address of result vector's data object.
       * @param  num_elems            Number of vectors's elements.
       */
      static inline
      const oclError &
      ocl_operator_equal              ( oclDataObject * const      arg1,
                                            elem_type            scalar,
                                        oclDataObject * const    result,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_equal (scalar)", op_v_level);
        ocl_basic_operator_kernel_21 ("scalar_equal", arg1, scalar, result, num_elems);
      
      }


      /**
       * @brief                       Elementwise equality of vectors.
       *
       * @param  arg1                 Address of first vector's data object.
       * @param  arg2                 Address of second vector's data object.
       * @param  result               Address of result vector's data object.
       * @param  num_elems            Number of vectors's elements.
       */
      static inline
      const oclError &
      ocl_operator_equal              ( oclDataObject * const      arg1,
                                        oclDataObject * const      arg2,
                                        oclDataObject * const    result,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_equal (vector)", op_v_level);
        ocl_basic_operator_kernel_3 ("vector_equal", arg1, arg2, result, num_elems);
      
      }
      
      
      /**
       * @brief                       Elementwise assignment of scalar to vector.
       *
       * @param  arg1                 Address of vector's data object.
       * @param  scalar               Scalar (to be assigned).
       * @param  num_elems            Number of vector's elements.
       */
      static inline
      const oclError &
      ocl_operator_assign             ( oclDataObject * const      arg1,
                                            elem_type            scalar,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_assign", op_v_level);
        ocl_basic_operator_kernel_11 ("assign", arg1, scalar, num_elems);
      
      } 
          
    
      //@}
      
      
      /**
       * @brief                       Deep copy of arg1 to arg2.
       *
       * @param  dest                 Address of destination.
       * @param  src                  Address of source.
       * @param  num_elems            Number of elements.
       */
      static inline
      const oclError &
      ocl_operator_copy               ( oclDataObject * const      dest,
                                        oclDataObject * const       src,
                                                  int         num_elems )
      {
      
        print_optional ("oclTraits <", trait :: print_elem_type (), "> :: ocl_operator_dec", op_v_level);
        ocl_basic_operator_kernel_2 ("copy_buffer", dest, src, num_elems);
      
      }
    
    
    
  };
  
  
  
  
# endif // __OCL_TRAITS_HPP__
