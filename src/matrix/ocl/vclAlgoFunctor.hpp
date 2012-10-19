# ifndef __VCL_ALGO_FUNCTOR_HPP__



  /************
   ** makros **
   ************/
  # define __VCL_ALGO_FUNCTOR_HPP__



  /**************************
   ** forward declarations **
   **************************/
  class oclViennaClObject;



  /***************************
   ** class: vclAlgoFunctor **
   **   (interface)         **
   **   (still containing   **
   **    dummy functions)   **
   ***************************/
  class vclAlgoFunctor
  {

  
    public:
  
      virtual
      void
      operator()              ()
      const
      {
        std::cout << "not implemented" << std::endl;
      }
    

    protected:
    
      // constructor
      vclAlgoFunctor          (const oclViennaClObject * const p_vclObj)
                            : mp_vclObj (p_vclObj)
      {
        std::cout << "Ctor: \"vclAlgoFunctor\"" << std::endl;
        /* TODO */
      }
      
      // member variables
      const oclViennaClObject * const mp_vclObj;


  }; // class vclAlgoFunctor
  
  
  
  /****************************
   ** class: vclSubtractAlgo **
   **    (derived)           **
   ****************************/
  class vclSubtractAlgo : public vclAlgoFunctor
  {
  
  
    public:
    
      // constructor
      vclSubtractAlgo       (const oclViennaClObject * const p_vclObj)
                          : vclAlgoFunctor (p_vclObj)
      {

        std::cout << "Ctor: \"vclSubtractAlgo\"" << std::endl;

        /* TODO */

      }
      
      virtual
      void
      operator()       ()
      const
      {

        std::cout << "vclSubtractAlgo :: operator() !!! float !!!" << std::endl;

//        viennacl :: vector <float> vclAlgoFunctor :: mp_vclObj -> getVCLArg < viennacl :: vector <float> > (1);

        /* TODO */
        
      }
  
  
  }; // class vclSubtractAlgo
  
    
  
  /**********************
   ** global functions **
   **********************/
  static
  const vclAlgoFunctor * const
  get_algo_functor     (const string & algo_name,
                        const oclViennaClObject * const p_vclObj)
  {
    
    std::cout << " :: get_algo_functor" << std::endl;
    std::cout << "!! choosing subtraction !!" << std::endl;
    
    return new vclSubtractAlgo (p_vclObj);
    
    /* TODO */
    
  }
  
  
  
# endif // __VCL_ALGO_FUNCTOR_HPP__
