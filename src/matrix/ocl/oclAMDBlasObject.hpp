# ifndef __OCL_AMD_BLAS_OBJECT_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_AMD_BLAS_OBJECT_HPP__
  
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclFunctionObject.hpp"
  # include "oclSettings.hpp"



  /**************************
   ** forward declarations **
   **************************/

  // classes
  template <class T, class S>
  class oclAMDBlasFunctor;
  template <class T, class S>
  class oclAMDBlasObject;

  // functions
  template <class T, class S>
  static
  const oclAMDBlasFunctor <T, S> * const
  get_amdBlas_functor     (const   oclAMDBlasType        &        algo_name,
                           const oclAMDBlasObject <T, S> * const   p_amdObj);




  /******************************
   ** class: oclAMDBlasObject  **
   **   (derived)              **
   ******************************/
  template <class T, class S>
  class oclAMDBlasObject : public oclFunctionObject
  {
  
  
    
    public:

      
      /**
       * @name            constructors and destructors
       */
      //@{
      
      
      /**
       * @brief           constructor
       */
      oclAMDBlasObject    ( const oclAMDBlasType                           algo,
                                   oclDataObject  * const * const   pp_vec_args,
                                             int                   num_vec_args,
                                             int                    num_scalars = 0,
                                             int  * const             p_scalars = NULL )
                         : oclFunctionObject (pp_vec_args, num_vec_args),
                           m_num_scalars     (num_scalars),
                           mp_scalars        (p_scalars),
                           mp_algo_functor   (get_amdBlas_functor <T> (algo, this))
      {
      
        print_optional ("Ctor: \"oclAMDBlasObject\"", VERB_HIGH);

        /* TODO */

      }
      
      
      /**
       * @brief           virtual destructor
       */
      virtual
      ~oclAMDBlasObject  ();
      
      
      //@}
      
      
      /**
       * @brief           run amd blas algorithm
       */
      virtual
      void
      run                 (const LaunchInformation &);
      
      
      /**
       * @brief           retrieve profiling information
       */
      virtual
      const ProfilingInformation
      getProfilingInformation   (const int)
      const;
      
      
      /**
       * @brief           retrieve size at given position
       */
      int
      getScalarArg        (const int num)
      const;
      
      
      /**
       * @brief           retrieve memory object of function argument
       */
      cl_mem
      getMemObject        (const int num)
      const;
      
      oclDataObject *
      getArg (const int num)
      const
      {
        return mpp_args [num];
      }
      
      
      
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      
      const               int                m_num_scalars;
                          int        * const mp_scalars;
    
      const oclAMDBlasFunctor <T, S> * const mp_algo_functor;
    
  
  }; // class oclAMDBlasObject
  
  
  
  
  /**************************
   ** function definitions **
   **************************/



  /**
   * @brief           retrieve size at given position
   */
  template <class T, class S>
  int
  oclAMDBlasObject <T, S> ::
  getScalarArg          (const int num)
  const
  {
  
    print_optional ("oclAMDBlasObject :: getScalarArg ()", VERB_HIGH);

    if (num > m_num_scalars)
    {
      throw oclError ("Requested size argument number is out of range!", "oclAMDBlasObject :: getScalarArg");
    }
    else if (mp_scalars == NULL)
    {
      throw oclError ("No Sizes given!", "oclAMDBlasObject :: getScalarArg");
    }
    
    return mp_scalars [num];
  
  }
  
  
  
  /**
   * @brief               refer to class definition
   */
  template <class T, class S>
  cl_mem
  oclAMDBlasObject <T, S> ::
  getMemObject            ( const int num )
  const
  {
    
    print_optional ("oclAMDBlasObject :: getMemObject ()", VERB_HIGH);
    
    if (num >= m_num_args)
    {
      throw oclError ("Requested argument number is out of range!", "oclAMDBlasObject :: getMemObject");
    }
    
    return (* (mpp_args [num] -> getBuffer ())) ();
    
  }




  /******************
   ** include (II) **
   ******************/
  # include "oclAMDBlasFunctor.hpp"



  /**************************************
   ** function definitions (continued) **
   **************************************/
  
  
  
  /**
   * @brief               virtual destructor
   */
  template <class T, class S>
  oclAMDBlasObject <T, S> ::
  ~oclAMDBlasObject      ()
  {
    
    print_optional ("Dtor: \"oclAMDBlasObject\"", VERB_HIGH);

    delete mp_algo_functor;

    /* TODO */

  }



  /**
   * @brief               refer to class definition
   */
  template <class T, class S>
  void
  oclAMDBlasObject <T, S> ::
  run                     (const LaunchInformation & lc)
  {
    
    print_optional ("oclAMDBlasObject :: run ()", VERB_HIGH);

    // prepare kernel arguments (load to gpu)
    for (int i = 0; i < m_num_args; i++)
    {

      // prepare argument
      mpp_args [i] -> prepare ();

    }
    
    // execute functor
    (*mp_algo_functor) ();
    
    // get data
    for (int i = 0; i < m_num_args; i++)
    {
      mpp_args [i] -> finish ();
    }
    
  }

  
  template <class T, class S>
  const ProfilingInformation
  oclAMDBlasObject <T, S> ::
  getProfilingInformation   (const int)
  const
  {
    throw oclError ("Not yet implemented!", "oclAMDBlasObject :: getProfilingInformation");
    
    return {-1.0, -1.0, -1.0, -1.0};
  }
  
  
  
# endif // __OCL_AMD_BLAS_OBJECT_HPP__
