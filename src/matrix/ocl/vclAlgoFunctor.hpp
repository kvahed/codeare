# ifndef __VCL_ALGO_FUNCTOR_HPP__




  /************
   ** makros **
   ************/
  # define __VCL_ALGO_FUNCTOR_HPP__



  /**************
   ** includes **
   **************/
   
  // ocl
  # include "oclSettings.hpp"




  /**************************
   ** forward declarations **
   **************************/
  template <class T>
  class oclViennaClObject;




  /***************************
   ** class: vclAlgoFunctor **
   **   (interface)         **
   **   (still containing   **
   **    dummy functions)   **
   ***************************/
  template <class T>
  class vclAlgoFunctor
  {
  
  

    public:


      /**
       * @brief               overloaded operator (functor)
       */
      virtual
      void
      operator()              ()
      const
      {
      
        std::cout << "not implemented" << std::endl;
        
      }
    
    

    protected:
    
    
      /**
       * @name                constructors and destructors
       */
      //@{
      
      
      /**
       * @brief               constructor
       */
      vclAlgoFunctor          ( const oclViennaClObject <T> * const p_vclObj )
                             : mp_vclObj (p_vclObj)
      {
      
        print_optional ("Ctor: \"vclAlgoFunctor\"", VERB_HIGH);
        
        /* TODO */
        
      }
      
      
      //@}
      
      
      /**********************
       ** member variables **
       **********************/
      
      const oclViennaClObject <T> * const mp_vclObj;


      /**********************
       ** static variables **
       **********************/
      static
      const VerbosityLevel op_v_level = VERB_MIDDLE;



  }; // class vclAlgoFunctor
  
  

  
  /****************************
   ** class: vclSubtractAlgo **
   **    (derived)           **
   ****************************/
  template <class T>
  class vclSubtractAlgo : public vclAlgoFunctor <T>
  {

  
  
    public:

    
      /**
       * @name              constructors and destructors
       */
      //@{
      
      
      /**
       * @brief             constructor
       */
      vclSubtractAlgo       ( const oclViennaClObject <T> * const p_vclObj )
                           : vclAlgoFunctor <T> (p_vclObj)
      {

        print_optional ("Ctor: \"vclSubtractAlgo\"", VERB_HIGH);

        /* TODO */

      }
      
      
      //@}
      
      
      /**
       * @brief             refer to base class
       */
      virtual
      void
      operator()            ()
      const
      {

        print_optional ("vclSubtractAlgo <T> :: operator()", vclAlgoFunctor <T> :: op_v_level);

        /**
         * create ViennaCl arguments
         */
        viennacl :: vector <T> arg1 = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (0);
        viennacl :: vector <T> arg2 = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (1);
        viennacl :: vector <T> diff = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (2);
        
        diff = arg1 - arg2;
        
//        std::cout << " -------------------\n arg1: \n" << arg1 << "\n -------------------" << std::endl;
//        std::cout << " -------------------\n arg2: \n" << arg2 << "\n -------------------" << std::endl;
//        std::cout << " -------------------\n diff: \n" << diff << "\n -------------------" << std::endl;

        /* TODO */
        
      }
  
  
  
  }; // class vclSubtractAlgo
  
  
  
  
  /**********************
   ** global functions **
   **********************/
  template <class T>
  static
  const vclAlgoFunctor <T> * const
  get_algo_functor     (const       vclAlgoType     &           algo,
                        const oclViennaClObject <T> * const p_vclObj)
  {
        
    print_optional (" :: get_algo_functor", VERB_HIGH);

    /* choose requested algorithm */
    switch (algo)
    {
    
      case vclSUBTRACT:
        return new vclSubtractAlgo <T> (p_vclObj);

      default:
        throw "*!* Specified ViennaCL algorithm not available! *!*";
    
    }

  }
  
  
  
  
# endif // __VCL_ALGO_FUNCTOR_HPP__
