# ifndef __OCL_DATA_OBJECT_HPP__



  # define __OCL_DATA_OBJECT_HPP__



  /**************
   ** includes **
   **************/

  // ocl
  # include "oclObservableDataObject.hpp"
  # include "oclSettings.hpp"
  # include "oclConnection.hpp"



  /**************************
   ** class: oclDataObject **
   **  (declaration)       **
   **************************/
  class oclDataObject : public oclObservableDataObject
  {


    public:
  
      // pure virtual: prepare ()
      virtual
      oclError
      prepare             () = 0;

      template <class T>
      T
      getVCLObject ()
      const
      {
        /* TODO */
      }


    protected:
      
      // constructor
      oclDataObject       ()
                        : m_gpu_obj_id (id_counter ++)
      {
        std::cout << "Ctor: \"oclDataObject\"" << std::endl;
        
        oclConnection :: Instance () -> addDataObject (this, m_gpu_obj_id);
        /* TODO */
      }
      
      const oclObjectID   m_gpu_obj_id;

    
    private:
    
      static oclObjectID  id_counter;


  };


  oclObjectID
  oclDataObject ::
  id_counter              = 0;


  
# endif // __OCL_DATA_OBJECT_HPP__
