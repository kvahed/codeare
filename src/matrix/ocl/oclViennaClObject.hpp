# ifndef __OCL_VIENNA_CL_OBJECT_HPP__




  /************
   ** makros **
   ************/
  # define __OCL_VIENNA_CL_OBJECT_HPP__

  
  
  
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
  class vclAlgoFunctor;
  template <class T, class S>
  class oclViennaClObject;

  // functions
  template <class T, class S>
  static
  const vclAlgoFunctor <T, S> * const
  get_algo_functor     (const       vclAlgoType &            algo_name,
                        const oclViennaClObject <T, S> * const   p_vclObj);
  
  
  
  
  /******************************
   ** class: oclViennaClObject **
   **   (derived)              **
   ******************************/
  template <class T, class S>
  class oclViennaClObject : public oclFunctionObject
  {
  
  
    
    public:

      
      /**
       * @name            constructors and destructors
       */
      //@{
      
      
      /**
       * @brief           constructor
       */
      oclViennaClObject   ( const   vclAlgoType                           algo,
                                  oclDataObject  * const * const   pp_vec_args,
                                            int                   num_vec_args,
                                            int                    num_scalars = 0,
                                            int  * const             p_scalars = NULL )
                         : oclFunctionObject (pp_vec_args, num_vec_args),
                           m_num_scalars     (num_scalars),
                           mp_scalars        (p_scalars),
                           mp_algo_functor   (get_algo_functor <T> (algo, this))
      {
      
        print_optional ("Ctor: \"oclViennaClObject\"", VERB_HIGH);

        /* TODO */

      }
      
      
      /**
       * @brief           virtual destructor
       */
      virtual
      ~oclViennaClObject  ();
      
      
      //@}
      
      
      /**
       * @brief           run vcl algorithm
       */
      virtual
      void
      run                 ();
      
      
      /**
       * @brief           retrieve size at given position
       */
      int
      getScalarArg          (const int num)
      const;
      
      
      /**
       * @brief           retrieve data wrapped by a ViennaCl matrix
       */
      template <class U>
      viennacl :: vector <U>
      getVCLVector        (const int num )
      const;
      

      /**
       * @brief           retrieve data wrapped by a ViennaCl matrix
       */
      template <class U,
                typename R>
      viennacl :: matrix <U, R>
      getVCLMatrix        (const        int  num,
                                        int    m,
                                        int    n )
      const;
      
      
      
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      
      const            int                m_num_scalars;
                       int        * const mp_scalars;
    
      const vclAlgoFunctor <T, S> * const mp_algo_functor;
    
  
  }; // class oclViennaClObject
  
  
  
  
  /**************************
   ** function definitions **
   **************************/



  /**
   * @brief           retrieve size at given position
   */
  template <class T, class S>
  int
  oclViennaClObject <T, S> ::
  getScalarArg          (const int num)
  const
  {
  
    print_optional ("oclViennaClObject :: getScalarArg ()", VERB_HIGH);

    if (num > m_num_scalars)
    {
      throw oclError ("Requested size argument number is out of range!", "oclViennaClObject :: getScalarArg");
    }
    else if (mp_scalars == NULL)
    {
      throw oclError ("No Sizes given!", "oclViennaClObject :: getScalarArg");
    }
    
    return mp_scalars [num];
  
  }
  
  
  
  /**
   * @brief               refer to class definition
   */
  template <class T, class S>
  template <class U>
  viennacl :: vector <U>
  oclViennaClObject <T, S> ::
  getVCLVector            ( const int num )
  const
  {
    
    print_optional ("oclViennaClObject :: getVCLVector ()", VERB_HIGH);
    
    if (num >= m_num_args)
    {
      throw oclError ("Requested argument number is out of range!", "oclViennaClObject :: getVCLVector");
    }
    
    return mpp_args [num] -> getVCLVector <U> ();
    
  }




  /**
   * @brief               refer to class definition
   */
  template <class T, class S>
  template <class U,
            typename R>// = viennacl :: column_major>
  viennacl :: matrix <U, R>
  oclViennaClObject <T, S> ::
  getVCLMatrix            ( const int num,
                            const int   m,
                            const int   n )
  const
  {
    
    print_optional ("oclViennaClObject :: getVCLMatrix ()", VERB_HIGH);
    
    if (num >= m_num_args)
    {
      throw oclError ("Requested argument number is out of range!", "oclViennaClObject :: getVCLMatrix");
    }
    
    try
    {
    
      return mpp_args [num] -> getVCLMatrix <U, R> (m, n);
  
    }
    catch (oclError & err)
    {
      throw oclError (err, "oclViennaClObject <T> :: getVCLMatrix");
    }
    
  }




  /******************
   ** include (II) **
   ******************/
  # include "vclAlgoFunctor.hpp"



  /**************************************
   ** function definitions (continued) **
   **************************************/
  
  
  
  /**
   * @brief               virtual destructor
   */
  template <class T, class S>
  oclViennaClObject <T, S> ::
  ~oclViennaClObject      ()
  {
    
    print_optional ("Dtor: \"oclViennaClObject\"", VERB_HIGH);

    delete mp_algo_functor;

    /* TODO */

  }



  /**
   * @brief               refer to class definition
   */
  template <class T, class S>
  void
  oclViennaClObject <T, S> ::
  run                     ()
  {
    
    print_optional ("oclViennaClObject :: run ()", VERB_HIGH);

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
  
  
  
  
# endif // __OCL_VIENNA_CL_OBJECT_HPP__
