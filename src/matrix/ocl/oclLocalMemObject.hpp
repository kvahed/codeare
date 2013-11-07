# ifndef __OCL_LOCAL_MEM_OBJECT_HPP__



  /************
   ** makros **
   ************/
  #define __OCL_LOCAL_MEM_OBJECT_HPP__
  
  
  
  /**************
   ** includes **
   **************/
  
  // ocl
  # include "oclSettings.hpp"
  # include "oclDataWrapper.hpp"
  
  
  
  /******************************
   ** class: oclLocalMemObject **
   **   (derived)              **
   ******************************/
  template <class T>
  class oclLocalMemObject : public oclDataWrapper <T>
  {
  
    
    
    public:
      
      
      /**
       * @name                constructors and destructors
       */
      //@{
      
      
      /**
       * @brief               constructor
       *
       * @param               cpu_data  @see oclDataWrapper<T>
       * @param               size      @see ocldataWrapper<T>
       *
       */
      oclLocalMemObject       (const    int   &        num_elems)
                            : oclDataWrapper <T> (NULL, num_elems)
      {
      
        print_optional ("Ctor: \"oclLocalMemObject\"", VERB_HIGH);
        
      }

      
      /**
       * @brief               virtual destructor
       */
      virtual
      ~oclLocalMemObject       ()
      {
      
        print_optional ("Dtor: \"oclLocalMemObject\"", VERB_HIGH);

      }

      
      //@}

    
      /**
       * @brief               inherited (oclDataObject)
       */
      virtual
      double
      prepare                 ();
      
      
      /**
       * @brief               inherited (oclObservableDataObject)
       */
      virtual
      double
      finish                  ();
      
      
      /**
       * @brief               inherited (oclDataWrapper)
       */
      virtual
      double
      getData                 ();
      
      
      /**
       * @brief               gpu upload size
       */
      virtual
      size_t
      getUploadSize           ()
      const;

    
     
  }; // class oclLocalMemObject
  
  
  
  /**************************
   ** function definitions **
   **************************/
  
  

  /**
   * @brief                   gpu upload size
   */
  template <class T>
  size_t
  oclLocalMemObject <T> ::
  getUploadSize               ()
  const
  {
    
    return oclDataObject :: getSize ();
    
  }
  
  

  
  
  /**
   * @brief                   do nothing
   */
  template <class T>
  double
  oclLocalMemObject <T> ::
  prepare                     ()
  {

    print_optional ("oclLocalMemObject::prepare (%d)", oclDataObject :: getID (), VERB_MIDDLE);
  
    // set status: calculating (set available via finish ())
    oclDataObject :: setLocked ();

    return .0;
    
  }
  
  
  
  /**
   * @brief                   do nothing
   */
  template <class T>
  double
  oclLocalMemObject <T> ::
  finish                      ()
  {

    print_optional ("oclLocalMemObject::finish", VERB_HIGH);

    // update data state: available for use
    oclDataObject :: setUnlocked ();

    return .0;
    
  }
  
  
  
  /**
   * @brief                   do nothing but throw error!
   */
  template <class T>
  double
  oclLocalMemObject <T> ::
  getData                     ()
  {
  
    print_optional ("oclLocalMemObject::getData", VERB_HIGH);
    
    throw oclError (" ** No Data! **", "oclLocalMemObject :: getData");
    
  }
  
  
  
# endif // __OCL_LOCAL_MEM_OBJECT_HPP__
