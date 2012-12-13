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
  template <class T>
  class vclAlgoFunctor;
  template <class T>
  class oclViennaClObject;

  // functions
  template <class T>
  static
  const vclAlgoFunctor <T> * const
  get_algo_functor     (const       vclAlgoType &            algo_name,
                        const oclViennaClObject <T> * const   p_vclObj);
  
  
  
  
  /******************************
   ** class: oclViennaClObject **
   **   (derived)              **
   ******************************/
  template <class T>
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
      template <class S>
      viennacl :: vector <S>
      getVCLVector        (const int num )
      const;
      

      /**
       * @brief           retrieve data wrapped by a ViennaCl matrix
       */
      template <class S,
                typename R>
      viennacl :: matrix <S, R>
      getVCLMatrix        (const        int  num,
                                        int    m,
                                        int    n )
      const;
      
      
      
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      
      const            int             m_num_scalars;
                       int     * const mp_scalars;
    
      const vclAlgoFunctor <T> * const mp_algo_functor;
    
  
  }; // class oclViennaClObject
  
  
  
  
  /**************************
   ** function definitions **
   **************************/



  /**
   * @brief           retrieve size at given position
   */
  template <class T>
  int
  oclViennaClObject <T> ::
  getScalarArg          (const int num)
  const
  {
  
    print_optional ("oclViennaClObject :: getSizeArg ()", VERB_HIGH);

    if (num > m_num_scalars)
    {
      throw oclError ("Requested size argument number is out of range!", "oclViennaClObject :: getSizeArg");
    }
/*    else if (mp_scalars == NULL)
    {
      throw oclError ("No Sizes given!", "oclViennaClObject :: getSizeArg");
    }
*/    
    return mp_scalars [num];
  
  }
  
  
  
  /**
   * @brief               refer to class definition
   */
  template <class T>
  template <class S>
  viennacl :: vector <S>
  oclViennaClObject <T> ::
  getVCLVector            ( const int num )
  const
  {
    
    print_optional ("oclViennaClObject :: getVCLArg ()", VERB_HIGH);
    
    if (num >= m_num_args)
    {
      throw oclError ("Requested argument number is out of range!", "oclViennaClObject :: getVCLVector");
    }
    
    return mpp_args [num] -> getVCLVector <S> ();
    
  }




  /**
   * @brief               refer to class definition
   */
  template <class T>
  template <class S,
            typename R>// = viennacl :: column_major>
  viennacl :: matrix <S, R>
  oclViennaClObject <T> ::
  getVCLMatrix            ( const int num,
                            const int   m,
                            const int   n )
  const
  {
    
    print_optional ("oclViennaClObject :: getVCLMatrix ()", VERB_HIGH);
    
    if (num >= m_num_args)
    {
      throw oclError ("Requested argument number is out of range!", "oclViennaClObject :: getVCLVector");
    }
    
    try
    {
    
      return mpp_args [num] -> getVCLMatrix <S, R> (m, n);
  
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
  template <class T>
  oclViennaClObject <T> ::
  ~oclViennaClObject      ()
  {
    
    print_optional ("Dtor: \"oclViennaClObject\"", VERB_HIGH);

    delete mp_algo_functor;

    /* TODO */

  }



  /**
   * @brief               refer to class definition
   */
  template <class T>
  void
  oclViennaClObject <T> ::
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
