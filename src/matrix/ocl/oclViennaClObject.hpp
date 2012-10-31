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


  
  /**************************
   ** forward declarations **
   **************************/

  // classes
  class vclAlgoFunctor;
  class oclViennaClObject;

  // functions
  static
  const vclAlgoFunctor * const
  get_algo_functor     (const string & algo_name,
                        const oclViennaClObject * const p_vclObj);
  
  
  
  
  /******************************
   ** class: oclViennaClObject **
   **   (derived)              **
   ******************************/
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
      oclViennaClObject   (const string                         algo_name,
                                 oclDataObject  * const * const pp_args,
                                 int                            num_args)
                        : oclFunctionObject (pp_args, num_args),
                          mp_algo_functor   (get_algo_functor (algo_name, this))
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
      template <class T>
      viennacl :: vector <T>
      getVCLArg           (const int num)
      const;
      
      
      
    protected:
    
    
      /**********************
       ** member variables **
       **********************/
      const vclAlgoFunctor * const mp_algo_functor;
    
    
  
  }; // class oclViennaClObject
  
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  
  /**
   * @brief               refer to class definition
   */
  template <class T>
  viennacl :: vector <T>
  oclViennaClObject ::
  getVCLArg               (const int num)
  const
  {
    
    std::cout << "oclViennaClObject :: getVCLArg" << std::endl;
    
    return mpp_args [num] -> getVCLObject <T> ();
    
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
  oclViennaClObject ::
  ~oclViennaClObject      ()
  {
    
    std::cout << "Dtor: \"oclViennaClObject\"" << std::endl;

    delete mp_algo_functor;

    /* TODO */

  }



  /**
   * @brief               refer to class definition
   */
  void
  oclViennaClObject ::
  run                     ()
  {
    
    std::cout << "oclViennaClObject :: run!" << std::endl;
    
    // execute functor
    (*mp_algo_functor) ();
    
  }
  
  
  
  
# endif // __OCL_VIENNA_CL_OBJECT_HPP__
