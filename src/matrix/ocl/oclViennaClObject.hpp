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
  # include "vclAlgoFunctor.hpp"
  
  
  
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
       * @brief           destructor
       */
      virtual
      ~oclViennaClObject  ()
      {
      
        std::cout << "Dtor: \"oclViennaClObject\"" << std::endl;

        delete mp_algo_functor;

        /* TODO */

      }
      
      //@}
      
      virtual
      void
      run                 ();
      
      template <class T, template <class S = T> class V >
      V <T>
      getVCLArg           (const int num)
      const;
      
      
    protected:
    
      const vclAlgoFunctor * const mp_algo_functor;
    
  
  }; // class oclViennaClObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  void
  oclViennaClObject ::
  run                     ()
  {
    
    std::cout << "oclViennaClObject :: run!" << std::endl;
    
    (*mp_algo_functor) ();
    
  }
  
  
  
  template <class T, template <class S = T> class V >
  V <T>
  oclViennaClObject ::
  getVCLArg               (const int num)
  const
  {
    
    std::cout << "oclViennaClObject :: getVCLArg" << std::endl;
    
    return mpp_args [num] -> getVCLObject < T, V<T> > ();
    
  }
  
  
  
# endif // __OCL_VIENNA_CL_OBJECT_HPP__
