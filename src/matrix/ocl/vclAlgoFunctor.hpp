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
  class vclAddSubAlgo : public vclAlgoFunctor <T>
  {

  
  
    public:

    
      /**
       * @name              constructors and destructors
       */
      //@{
      
      
      /**
       * @brief             constructor
       */
      vclAddSubAlgo       ( const oclViennaClObject <T> * const p_vclObj,
                            const               int                 sign )
                           : vclAlgoFunctor <T> (p_vclObj),
                             m_sign                 (sign)
      {

        print_optional ("Ctor: \"vclAddSubAlgo\"", VERB_HIGH);

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

        print_optional ("vclAddSubAlgo <T> :: operator()", vclAlgoFunctor <T> :: op_v_level);

        /**
         * create ViennaCl arguments
         */
        viennacl :: vector <T>   arg1 = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (0);
        viennacl :: vector <T>   arg2 = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (1);
        viennacl :: vector <T> result = vclAlgoFunctor <T> :: mp_vclObj -> template getVCLArg <T> (2);
        
        /* perform calculation */
        if (m_sign < 0)
          result = arg1 - arg2;
        else
          result = arg1 + arg2;
        
      }
      
      
    private:
    
      int m_sign;
  
  
  
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
        return new vclAddSubAlgo <T> (p_vclObj, -1);

      default:
        throw "*!* Specified ViennaCL algorithm not available! *!*";
    
    }

  }
  
  
  
  
# endif // __VCL_ALGO_FUNCTOR_HPP__
