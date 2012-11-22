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
  get_algo_functor     (const       vclAlgoType &       algo_name,
                        const oclViennaClObject <T> * const  p_vclObj);
  
  
  
  
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
      oclViennaClObject   (const vclAlgoType                        algo,
                                 oclDataObject  * const * const  pp_args,
                                 int                            num_args)
                        : oclFunctionObject (pp_args, num_args),
                          mp_algo_functor   (get_algo_functor <T> (algo, this))
      {
      
        std::cout << "Ctor: \"oclViennaClObject\"" << std::endl;

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
       * @brief           retrieve data wrapped by a ViennaCl object
       */
      template <class S>
      viennacl :: vector <S>
      getVCLArg           (const int num)
      const;
      
      
      
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      const vclAlgoFunctor <T> * const mp_algo_functor;
    
    
  
  }; // class oclViennaClObject
  
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  
  /**
   * @brief               refer to class definition
   */
  template <class T>
  template <class S>
  viennacl :: vector <S>
  oclViennaClObject <T> ::
  getVCLArg               (const int num)
  const
  {
    
    std::cout << "oclViennaClObject :: getVCLArg" << std::endl;
    
    return mpp_args [num] -> getVCLObject <S> ();
    
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
    
    std::cout << "Dtor: \"oclViennaClObject\"" << std::endl;

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
    
    std::cout << "oclViennaClObject :: run!" << std::endl;

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
