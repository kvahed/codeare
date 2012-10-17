# ifndef __OCL_OBSERVABLE_DATA_OBJECT_HPP__



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

      // pure virtual: finish ()
      virtual
      oclError
      finish () = 0;


  };



# endif // __OCL_OBSERVABLE_DATA_OBJECT_HPP__
