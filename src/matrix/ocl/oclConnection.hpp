# ifndef __OCL_CONNECTION_HPP__



  /************
   ** makros **
   ************/
  # define __OCL_CONNECTION_HPP__ 



  /**************
   ** includes **
   **************/

  // C++ std. headers
  # include <map>



  // ocl
  # include "oclSettings.hpp"



  /**************************
   ** forward declarations **
   **************************/
  class oclDataObject;
  class oclFunctionObject;
  class vclAlgoFunctor;



  /**************************
   ** class: oclConnection **
   **  (singleton)         **
   **************************/
  class oclConnection
  {

  
    public:

      /***********************
       ** enumeration types **
       ***********************/
      enum KernelType {KERNEL, VCL};
      enum SyncType   {SYNC, ASYNC};

      /**********************
       ** type definitions **
       **********************/
      

      /***************
       ** functions **
       ***************/
      static
      oclConnection *
      Instance              ();

      void
      Initialise            ();
    
      void
      Finalise              ();
      
      const oclFunctionObject * const
      makeFunctionObject    (const string & kernel_name,
                             oclDataObject * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type);
    
      const oclFunctionObject * const
      makeFunctionObject    (const vclAlgoFunctor & algo,
                             oclDataObject * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type);
                          
      void
      addDataObject         (oclDataObject * const ocl_obj, const oclObjectID & obj_id);
                          
      void
      run                   (oclFunctionObject * const func_obj) const;
    
    
    private:
    
      /**********************
       ** member variables **
       **********************/
      static oclConnection                    * mp_inst;
      std::map <oclObjectID, oclDataObject *>   m_current_ocl_objects;
  
      /******************
       ** constructors **
       ******************/
      oclConnection         () {/*TODO*/};
      oclConnection         (oclConnection &) { /*TODO*/};
  

  }; // class oclConnection
  
  

  /*******************
   ** includes (II) **
   *******************/
   
  // ocl
  # include "oclFunctionObject.hpp"
  # include "vclAlgoFunctor.hpp"
  # include "oclDataObject.hpp"
  
  

  /**************************************
   ** initialization of static members **
   **************************************/  
  oclConnection *
  oclConnection ::
  mp_inst                   = NULL;
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  
  oclConnection *
  oclConnection ::
  Instance ()
  {
  
    if (mp_inst == NULL)
      mp_inst = new oclConnection ();
      
    return mp_inst;
    
  }
  
  
  void
  oclConnection ::
  addDataObject (oclDataObject * const ocl_obj, const oclObjectID & obj_id)
  {
    std::cout << "addDataObject" << std::endl;
    m_current_ocl_objects.insert (std::pair <oclObjectID, oclDataObject *> (obj_id, ocl_obj));
  }



# endif // __OCL_CONNECTION_HPP__
