# ifndef __OCL_DATA_OBJECT_HPP__



  /************
   ** makros **
   ************/
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

  
      /**
       * @brief             pure virtual: prepare ()
       *
       *                    -- prepare data object for use on gpu --
       */
      virtual
      oclError &
      prepare               () = 0;


      /**
       * @brief             -- return ViennaCl representation of data object --
       */
      template <class S>
      viennacl :: vector <S>
      getVCLObject          ();


      /**
       * @name              getters (public)
       */
      //@{


      /**
       * @brief             return object's ID
       */
      inline
      oclObjectID
      getID                 () const;


      /**
       * @brief             return object's size
       */
      inline
      size_t
      getSize               () const;
      
      
      /**
       * @brief             return sync status (GPU <-> CPU memory)
       */
      inline
      bool
      getSyncState          () const;
      
      
      /**
       * @brief             is data available on GPU ?
       */
      inline
      bool
      getMemState           () const;
      
      
      /**
       * @brief             getter for buffer pointer
       */
      inline
      cl::Buffer *
      getBuffer             () const;
      
      
      //@}
      
      
      /**
       * @name              setters (public)
       */
      //@{
      
      
      /**
       * @brief             tell object that data (CPU) was modified
       */
      inline
      void
      setCPUModified        ();
      
      
      /**
       * @brief             set a buffer object (to be used by oclConnection)
       */
      inline
      void
      setBuffer             (cl::Buffer * p_gpu_buffer);
            
      
      //@}
      
      
      /**
       * @name              constructors and destructors
       */
      //@{

      
      /**
       * @brief             default constructor
       *
       * @param             size - ... in bytes
       */
      oclDataObject         (const size_t size)
                          : m_gpu_obj_id  (id_counter ++),      // init object ID
                            mp_gpu_buffer (NULL),               // init buffer (pointer: NULL)
                            m_size        (size),               // init size
                            m_on_gpu      (false),              // data loaded to GPU on demand
                            mp_modified   ({true, false}),      // ... so data aren't sync
                            m_lock        (false)               // ... and there're no calculations on GPU
      {
      
        std::cout << "Ctor: \"oclDataObject\" (" << m_gpu_obj_id << ")" << std::endl;
        
        // register data object at oclConnection
        oclConnection :: Instance () -> addDataObject (this);
        
      }
      

      /**
       * @brief             virtual destructor
       */
      virtual
      ~oclDataObject        ()
      {
      
        std::cout << "Dtor: \"oclDataObject\" (" << m_gpu_obj_id << ")" << std::endl;
        
        // unregister data object at oclConnection
        oclConnection :: Instance () -> removeDataObject (this);
        
        // delete buffer object
        delete mp_gpu_buffer;
        
      }

      
      //@}
      
      

    protected:


      /**
       * @brief             synchronize data on GPU with CPU
       */
      virtual
      oclError &
      loadToGPU             () = 0;


      /**
       * @brief             synchronize data on CPU with GPU
       */
      virtual
      oclError &
      loadToCPU             () = 0;


      /**
       * @name              setters (protected)
       */
      //@{
      
      
      /**
       * @brief             tell object that data is loaded to GPU
       */
      inline
      void
      setLoaded             ();
      
      
      /**
       * @brief             block data
       */
      inline
      void
      setLocked             ();
      
      
      /**
       * @brief             unblock data
       */
      inline
      void
      setUnlocked           ();
      
      
      /**
       * @brief             notify modification of data in GPU memory
       */
      inline
      void
      setGPUModified        ();
      
      
      /**
       * @brief             notify synchronicity of CPU and GPU data
       */
      inline
      void
      setSync               ();
      
      
      //@}
      
      
      /**
       * @name              getters (protected)
       */
      //@{
      
      
      /**
       * @brief             getter for m_busy
       */
      inline
      bool
      getLockState          () const;
      
      
      //@}
      
            
      /**********************
       ** member variables **
       **   (protected)    **
       **********************/
      const oclObjectID     m_gpu_obj_id;     // local id of particular object
             cl::Buffer   * mp_gpu_buffer;    // Buffer representation of data object
      
                 size_t     m_size;           // size of data (in bytes)
                 
                   bool     m_on_gpu;         // determines wether data is available on GPU
                   bool     mp_modified [2];  // specifies, if data on GPU and CPU are the same
                   bool     m_lock;           // determines wether data is used on GPU right now
    


    private:

    
      static oclObjectID    id_counter;       // global counter to produce unique IDs

      /****************
       ** enum types **
       ****************/
      
      // indices of m_modified
      enum {CPU, GPU};
      


  }; // class oclDataObject



  /**
   *  initialise static class member
   */
  oclObjectID
  oclDataObject ::
  id_counter              = 0;



  /**************************
   ** function definitions **
   **************************/



  /**
   *                      -- refer to class definition --
   */
  template <class S>
  viennacl :: vector <S>
  oclDataObject ::
  getVCLObject            ()
  {
  
    std::cout << "oclDataObject :: getVCLObject (!!! not yet implemented !!!)" << std::endl;
  
//    loadToGPU ();
  
//    std::cout << " ** after loadToGPU ()" << std::endl;
  
//    return viennacl :: vector <S> ((*mp_gpu_buffer)(), m_size);
  
    return viennacl :: vector <S> (1);
  
    /* TODO */
    
  }

  
  
        /*************
    --   ** getters **   --
         *************/
  
  
  /**
   * @brief               getter for member buffer
   */
  inline
  cl::Buffer *
  oclDataObject ::
  getBuffer               ()
  const
  {
  
    return mp_gpu_buffer;
  
  }
  
  
  /**
   * @brief               getter for object ID
   */
  inline
  oclObjectID
  oclDataObject ::
  getID                   ()
  const
  {
    
    return m_gpu_obj_id;
    
  }
  
  
  
  /**
   * @brief               getter of object size
   */
  inline
  size_t
  oclDataObject ::
  getSize                 ()
  const
  {
    
    return m_size;
    
  }
  
  
  
  /**
   * @brief               getter for memory state
   */
  inline
  bool
  oclDataObject ::
  getMemState             ()
  const
  {
  
    return m_on_gpu;
  
  }
  
  
  
  /**
   * @brief               getter for sync state
   */
  inline
  bool
  oclDataObject ::
  getSyncState            ()
  const
  {
    
    return     mp_modified [CPU]
            && mp_modified [GPU];
    
  }
  
 
  
  /**
   * @brief               getter for m_lock
   */
  inline
  bool
  oclDataObject ::
  getLockState            ()
  const
  {
  
    return m_lock;
    
  }
  
  
  
        /*************
    --   ** setters **   --
         *************/
  
  
  
  /**
   * @brief               set CPU data modified
   */
  inline
  void
  oclDataObject ::
  setCPUModified          ()
  {
    
    mp_modified [CPU] = true;
    
  }
  
  
  
  /**
   * @brief               set GPU data modified
   */
  inline
  void
  oclDataObject ::
  setGPUModified          ()
  {
    
    mp_modified [GPU] = true;
    
  }
  
  
  
  /**
   * @brief               set m_on_gpu
   */
  inline
  void
  oclDataObject ::
  setLoaded               ()
  {
  
    m_on_gpu = true;
  
  }
  
  
  
  /**
   * @brief               lock data (no object internal access)
   */
  inline
  void
  oclDataObject ::
  setLocked               ()
  {
  
    m_lock = true;
  
  }
  
  
  
  /**
   * @brief               unlock data
   */
  inline
  void
  oclDataObject ::
  setUnlocked             ()
  {
  
    m_lock = false;
  
  }
  
  
  
  /**
   * @brief               notify synchronicity
   */
  inline
  void
  oclDataObject ::
  setSync                 ()
  {
  
      mp_modified [CPU]
    = mp_modified [GPU]
    = false;
  
  }
  
  
  
  /**
   * @brief                set buffer pointer
   */
  inline
  void
  oclDataObject ::
  setBuffer             (cl::Buffer * p_gpu_buffer)
  {
  
    if (! getMemState ())
    {
      mp_gpu_buffer = p_gpu_buffer;
    }
    else
    {
      std::cout << " *!* Caution: Buffer already exists! *!*" << std::endl;
      throw int (-1);
    }
  
  }
   
   
  
# endif // __OCL_DATA_OBJECT_HPP__
