# ifndef __OCL_OBSERVABLE_DATA_OBJECT_HPP__



  /************
   ** makros **
   ************/
  # define __OCL_OBSERVABLE_DATA_OBJECT_HPP__
  
  

  /**************
   ** includes **
   **************/

  // ocl - general
  # include "oclError.hpp"



  /************************************
   ** class: oclObservableDataObject **
   **  (declaration)                 **
   ************************************/
  class oclObservableDataObject
  {


    public:

      /**
       * @brief            pure virtual: finish ()
       */
      virtual
      oclError &
      finish               () = 0;


  }; // class oclObservableDataObject



# endif // __OCL_OBSERVABLE_DATA_OBJECT_HPP__
