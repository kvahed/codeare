# ifndef __OCL_DATA_OBJECT_HPP__



  # define __OCL_DATA_OBJECT_HPP__



  /**************
   ** includes **
   **************/

  // ocl
  # include "oclObservableDataObject.hpp"
  # include "oclSettings.hpp"
  # include "oclConnection.hpp"

  // ViennaCL
  # include "/usr/include/viennacl/vector.hpp"



  /**************************
   ** class: oclDataObject **
   **  (declaration)       **
   **************************/
  class oclDataObject : public oclObservableDataObject
  {


    public:
  
      // pure virtual: prepare ()
      virtual
      oclError &
      prepare             (const int num) = 0;

      template <class T, template <class S = T> class V>
      V <T>
      getVCLObject        ()
      const;

      // constructor
      oclDataObject       ()
                        : m_gpu_obj_id (id_counter ++)
      {
        std::cout << "Ctor: \"oclDataObject\"" << std::endl;
        
        oclConnection :: Instance () -> addDataObject (this, m_gpu_obj_id);
        /* TODO */
      }
      
      // destructor
      virtual
      ~oclDataObject      ()
      {
        std::cout << "Dtor: \"oclDataObject\"" << std::endl;
        /* TODO */
      }


    protected:
      
      const oclObjectID   m_gpu_obj_id;

    
    private:
    
      static oclObjectID  id_counter;


  }; // class oclDataObject



  oclObjectID
  oclDataObject ::
  id_counter              = 0;



  /**************************
   ** function definitions **
   **************************/
  template <class T, template <class S = T> class V>
  V <T>
  oclDataObject ::
  getVCLObject            ()
  const
  {
  
    std::cout << "oclDataObject :: getVCLObject" << std::endl;
  
    return V <T> ();
  
    /* TODO */
    
  }


  
# endif // __OCL_DATA_OBJECT_HPP__
