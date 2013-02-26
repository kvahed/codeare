# ifndef __OCL_CONNECTION_HPP__



  /************
   ** makros **
   ************/
  # define __OCL_CONNECTION_HPP__ 
  # define __CL_ENABLE_EXCEPTIONS



  /**************
   ** includes **
   **************/

  // OpenCL
  # include <CL/cl.hpp>

  // C++ std. headers
  # include <map>
  # include <list>

  // ocl
  # include "oclSettings.hpp"
  # include "oclError.hpp"



  /**************************
   ** forward declarations **
   **************************/
  class oclDataObject;
  class oclFunctionObject;
//  class vclAlgoFunctor;
  


inline
std::string
modify_kernel        ( std::string const &       source,
                       std::string const & fp_extension,
                       std::string const &  vector_type,
                       std::string const &  scalar_type,
                               int const &            n = -1 );



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
      enum KernelType {KERNEL, VCL, AMD};
      enum SyncType   {SYNC, ASYNC};

      /**********************
       ** type definitions **
       **********************/
      typedef cl::Platform      clPlatform;
      typedef cl::Device        clDevice;
      typedef cl::Context       clContext;
      typedef cl::Program       clProgram;
      typedef cl::CommandQueue  clCommandQueue;
      typedef cl::Kernel        clKernel;
      typedef cl::Buffer        clBuffer;
      typedef std::vector <clPlatform>      clPlatforms;
      typedef std::vector <clDevice>        clDevices;
      typedef std::vector <clCommandQueue>  clCommandQueues;
      typedef std::vector <clKernel>        clKernels;
      typedef std::vector <clBuffer *>      clBuffers;
      

      /***************
       ** functions **
       ***************/
      inline
      static
      oclConnection *
      Instance              ();

      void
      Initialise            ();
    
      void
      Finalise              ();
      
      template <class T, class S>
      oclFunctionObject * const
      makeFunctionObject    (const string & kernel_name,
                             oclDataObject * const * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type);
                        
      template <class T, class S>     
      oclFunctionObject * const
      makeFunctionObject    (const vclAlgoType & algo_name,
                             oclDataObject * const * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type,
                             const int num_sclars = 0,
                             int * scalars = NULL);

      template <class T, class S>     
      oclFunctionObject * const
      makeFunctionObject    (const oclAMDBlasType & algo_name,
                             oclDataObject * const * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type,
                             const int num_sclars = 0,
                             int * scalars = NULL);
                             
 
      void
      addDataObject         (oclDataObject * const ocl_obj);
      
      void
      removeDataObject      (const oclDataObject * const ocl_obj);
      void
      releaseBuffer         (const oclDataObject * const ocl_obj);
                          
      void
      run                   (oclFunctionObject * const func_obj) const;

      // getters for OpenCL objects
      const clPlatform      getPlatform       () const { return m_plat; }
      const clDevices       getDevices        () const { return m_devs; }
      const clContext       getContext        () const { return m_cont; }
      const clCommandQueues getCommandQueues  () const { return m_comqs; }
      template <class T, class S>
      const clProgram       getProgram        () const { return ocl_precision_trait <T, S> :: getProgram (this); }
      
      
      /****************************************
       ** functions for OpenCL functionality **
       ****************************************/
      
      /**
       * @brief             load cpu data to GPU and return corresponding buffer
       *                      
       */
      template <class T>
      void
      loadToGPU             (       T   * const cpu_arg,
                             ::size_t           size,
                             clBuffer   * const buffer);
      
      
      /**
       * @brief             create buffer for given object id
       */
      template <class T>
      void
      createBuffer          (                T   * const cpu_arg,
                             const    ::size_t           size,
                             const oclObjectID           obj_id);
      
      
      /**
       * @brief             load data from given buffer to given cpu pointer
       *
       * @param             size - ... in bytes
       */
      template <class T>
      void
      loadToCPU             (const clBuffer * const buffer,
                                          T * const cpu_arg,
                             const ::size_t         size)
      {
      
        print_optional ("oclConnection :: loadToCPU ()", m_verbose);
        
        try
        {
        
          // create new buffer
          clBuffer * p_tmp_buffer = new clBuffer(m_cont, CL_MEM_READ_WRITE, size, cpu_arg, &m_error);
          
          print_optional (" -- size: %d Bytes", size, VERB_LOW); //m_verbose);
          
          // read data from given buffer
          m_error = m_comqs [0] . enqueueReadBuffer (*buffer, CL_TRUE, 0, size, cpu_arg, NULL, NULL);
        
        }
        catch (cl::Error cle)
        {
          
          throw oclError (cle.err(), "oclConnection :: loadToCPU");
          
        }
      
      }
      
  
      // activate kernel (stays activated until another kernel is activated
//      template <class T>
      int
      activateKernel        (const std::string kernelname);
    
      // run kernel with given dimensions
      int
      runKernel             (const cl::NDRange & global_dims,
                             const cl::NDRange & local_dims);
  
      
      
      /**
       * @brief               register kernel argument
       */
      cl_int
      setKernelArg            (      int              num,
                               const oclDataObject *  obj_id);
        
    
      void
      print_ocl_objects      ();
    

      // set central vars for particular precision mode for type T
      template <class T, class S>
      void
      activate              ( );
      
      
      // for amd blas
      cl_command_queue
      getCommandQueue       ( )
      const
      {
      
        return (m_comqs [0]) ();
      
      }


    
    /*************
     ** private **
     *************/
    private:


      /**********************
       ** member variables **
       **********************/
      static oclConnection                                * mp_inst;
      std::map <oclObjectID, oclDataObject * const>         m_current_ocl_objects;
      std::list <oclObjectID>                               m_loaded_ocl_objects;
      std::map <clBuffer *, int>                            m_current_buffers;
      
      
      /********************
       ** OpenCL members **
       ********************/
      clPlatform               m_plat;       // available platforms
      clDevices                m_devs;       // vector of available devices on m_plat
      clContext                m_cont;       // context used associated with m_devs
      clCommandQueues          m_comqs;      // each command queue is associated with corresponding device from m_devs
      clProgram              * mp_prog,
                               m_prog_ff,     // program containing whole source code (single precision)
                               m_prog_dd,     // program containing whole source code (double precision)
                               m_prog_cfcf,
                               m_prog_cff,
                               m_prog_cdcd,
                               m_prog_cdd;
      clKernels              * mp_kernels,
                               m_kernels_ff,  // kernels available in m_prog_f
                               m_kernels_dd,  // kernels available in m_prog_d
                               m_kernels_cfcf,
                               m_kernels_cff,
                               m_kernels_cdcd,
                               m_kernels_cdd;
//      std::vector<clBuffers> * mp_buffers,
//                               m_buffers_f,  // vectors of buffers for each kernels arguments (single prec)
//                               m_buffers_d;  // vectors of buffers for each kernels arguments (double prec)
      cl_int                   m_error;      // error variable
      VerbosityLevel           m_verbose;    // verbosity level
      clKernel               * mp_actKernel;
      int                      num_kernel;

      
      /*********************
       ** precision trait **
       **   (base struct) **
       *********************/
      template <class T, class S>
      struct ocl_precision_trait
      { };


      template <class T, class S> //, class trait = read_kernel_source_trait <T> >
      static const char * ReadSource   (std::string, int * size);
             int          BuildProgram ();
             const char * errorString (cl_int e);
             
      template <class T, class S>
      static
      oclError &
      init_program_kernels    (oclConnection * const);

  
      /******************
       ** constructors **
       ******************/
      oclConnection         (const char *,
                             const char *,
                             const char *,
                             const char *,
                             const char *,
                             const char *,
                             cl_device_type device_type = CL_DEVICE_TYPE_GPU,
                             VerbosityLevel verbose = VERB_NONE);
      
      oclConnection         (oclConnection &) { /*TODO*/};
    
      

  }; // class oclConnection
  
  
  
  /** ********************************************* **
   * @name            precision traits for kernels  **
   *                  in oclConnection              **
   ** ********************************************* **/
  //@{
  
  
  /**
   * forbidden type combination
   */
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <cxfl, double>
  { /* -- */ };
  
  
  /***********************************
   ** precision trait (cxfl, float) **
   **   (derived)                   **
   ***********************************/
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <cxfl, float>
  {
 
 
    public:
  
      typedef cxfl elem_type;
       
      static inline
      const char *
      getTypeString ()
      { return "<<cxfl, float>>"; }

      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/mixed_A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/mixed_AB_kernels.cl"));
        return filenames;
      }
          
      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: disable"),
                                                            std::string ("float2"),
                                                            std::string ("float"),
                                                                      8 ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_cff); }
            
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_cff); }

  };


  /***********************************
   ** precision trait (cxfl, cxfl)  **
   **   (derived)                   **
   ***********************************/
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <cxfl, cxfl>
  {
 
 
    public:
  
      typedef cxfl elem_type;
      
      static inline
      const char *
      getTypeString ()
      { return "<<cxfl, cxfl>>"; }

      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/complex_A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/complex_AB_kernels.cl"));
        return filenames;
      }
          
      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: disable"),
                                                            std::string ("float2"),
                                                            std::string ("float2") ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_cfcf); }
            
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_cfcf); }

  };      


  /************************************
   ** precision trait (float, float) **
   **   (derived)                    **
   ************************************/  
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <float, float>
  {

    public:

      typedef float elem_type;
        
      static inline
      const char *
      getTypeString ()
      { return "<<float, float>>"; }
          
      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/AB_kernels.cl"));
        return filenames;
      }

      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: disable"),
                                                            std::string ("float"),
                                                            std::string ("float"),
                                                                      8 ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_ff); }
         
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_ff); }

  };


  /**************************************
   ** precision trait (double, double) **
   **   (derived)                      **
   **************************************/
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <double, double>
  {
  
    public:

      typedef double elem_type;
       
      static inline
      const char *
      getTypeString ()
      { return "<<double, double>>"; }
          
      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/AB_kernels.cl"));
        return filenames;
      }

      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: enable"),
                                                            std::string ("double"),
                                                            std::string ("double"),
                                                                      8 ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_dd); }
          
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_dd); }
          
  };


  /************************************
   ** precision trait (cxdb, double) **
   **   (derived)                    **
   ************************************/
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <cxdb, double>
  {

    public:
 

      typedef double elem_type;
        
      static inline
      const char *
      getTypeString ()
      { return "<<cxdb, double>>"; }
          
      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/mixed_A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/mixed_AB_kernels.cl"));
        return filenames;
      }

      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: enable"),
                                                            std::string ("double2"),
                                                            std::string ("double"),
                                                                      8 ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_cdd); }
         
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_cdd); }

  };


  /***********************************
   ** precision trait (cxdb, cxdb)  **
   **   (derived)                   **
   ***********************************/
  template <>
  template <>
  struct oclConnection :: ocl_precision_trait <cxdb, cxdb>
  {
   
    public:
  
      typedef  float elem_type_A;
      typedef double elem_type_B;

      static inline
      const char *
      getTypeString ()
      { return "<<cxdb, cxdb>>"; }
  
      static inline
      const std::vector <std::string>
      getFilenames  ()
      {
        std::vector <std::string> filenames;
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/complex_A_kernels.cl"));
        filenames.push_back (std::string ("./src/matrix/ocl/kernels/complex_AB_kernels.cl"));
        return filenames;
      }

      static inline
      const char *
      modify_source (void * buf, int * size)
      {
        std::string * tmp_str = new string ( modify_kernel (std::string ((const char *) buf),
                                                            std::string ("cl_khr_fp64: enable"),
                                                            std::string ("double2"),
                                                            std::string ("double2") ) );
        *size = tmp_str -> size () * sizeof (char);
        return tmp_str->c_str ();
      }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_cdcd); }
        
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_cdcd); }
          
  };


  /***********************************
   ** precision trait (bool, S)     **
   **   (derived)                   **
   ***********************************/
  template <>
  template <class S>
  struct oclConnection :: ocl_precision_trait <bool, S>
  {

    public:

      typedef bool elem_type;

      static inline
      const std::vector <const char *>
      getFilenames  ()
      { return std::vector <const char *> (0); }

      static inline
      const char *
      modify_source (void * buf, int * size)
      { return (const char *) buf; }
          
      static inline
      clProgram *
      getProgram    (oclConnection * const con)
      { return & (con -> m_prog_ff); }
        
      static inline
      clKernels *
      getKernels    (oclConnection * const con)
      { return & (con -> m_kernels_ff); }

  };
  
  
  //@}
  
  

  /*******************
   ** includes (II) **
   *******************/
   
  // ocl
  # include "oclFunctionObject.hpp"
  # include "oclDataObject.hpp"
  # include "oclKernelObject.hpp"
  # include "oclViennaClObject.hpp"
  # include "oclAsyncKernelObject.hpp"
  # include "oclAMDBlasObject.hpp"
  
  

  /**************************************
   ** initialization of static members **
   **************************************/  
  oclConnection *
  oclConnection ::
  mp_inst                   = NULL;
  
  
  
  /**************************
   ** function definitions **
   **************************/

  # include "oclConnection.cpp"



  // set central vars for particular precision mode for type T (vector) and type S (scalar)
  template <class T, class S>
  void
  oclConnection ::
  activate              ( )
  {
    mp_prog    = ocl_precision_trait <T, S> :: getProgram (this);
    mp_kernels = ocl_precision_trait <T, S> :: getKernels (this);
  }

  


  void
  oclConnection :: 
  print_ocl_objects      ()
  {
      
    std::cout << " --------------------------------------------------- " << std::endl;
    std::cout << " oclDataObjects:" << std::endl;
    std::map <oclObjectID, oclDataObject * const> :: iterator it;
    for (it = m_current_ocl_objects.begin (); it != m_current_ocl_objects.end (); it++)
    {
        
      std::cout << " (" << it -> second -> getID () << ")" << std::endl;
        
    }
      
    std::cout << " --------------------------------------------------- " << std::endl;
      
  }


  
  /**
   * @brief             retrieve instance of oclConnection
   */
  inline
  oclConnection *
  oclConnection ::
  Instance              ()
  {
  
    if (mp_inst == NULL)
      mp_inst = new oclConnection ( "./src/matrix/ocl/kernels/A_kernels.cl",
                                    "./src/matrix/ocl/kernels/AB_kernels.cl",
                                    "./src/matrix/ocl/kernels/complex_A_kernels.cl",
                                    "./src/matrix/ocl/kernels/complex_AB_kernels.cl",
                                    "./src/matrix/ocl/kernels/mixed_A_kernels.cl",
                                    "./src/matrix/ocl/kernels/mixed_AB_kernels.cl",
                                    CL_DEVICE_TYPE_GPU,
                                    global_verbosity [OCL_CONNECTION] );
      
    return mp_inst;
    
  }
  
  
  
  /**
   * @brief               register kernel argument
   */
  cl_int
  oclConnection ::
  setKernelArg            (      int              num,
                           const oclDataObject *  p_arg)
  {
  
    // register kernel argument
    m_error = mp_actKernel -> setArg (num, * p_arg -> getBuffer ());
    
    print_optional (" -!-> setKernelArg: ", errorString (m_error), m_verbose);
  
  }
  
  
  
  /**
   * @brief             add oclDataObject to list of existing objects
   */
  void
  oclConnection ::
  addDataObject         (oclDataObject * const ocl_obj)
  {
  
    print_optional ("oclConnection :: addDataObject (...)", m_verbose);

    // insert new data object
    m_current_ocl_objects.insert (std::pair <oclObjectID, oclDataObject * const> (ocl_obj -> getID (), ocl_obj));
    
    // in case of existing buffer, update corresponding lists
    if (ocl_obj -> getMemState ())
    {
      
      print_optional (" *!* existing buffer (in addDObj) *!* ", m_verbose);
      
      // list of loaded (on GPU) oclObjects
      m_loaded_ocl_objects.push_back (ocl_obj -> getID ());
      
      // increase reference count for buffer
      std::map <clBuffer *, int> :: iterator tmp_it = m_current_buffers.find (ocl_obj -> getBuffer ());
      
      if (tmp_it == m_current_buffers.end ())
      {
      
        throw oclError ("Buffer not found", "oclConnection :: addDataObject");
        
      }

      // increase
      tmp_it -> second ++;
      
    }

  }
  
  
  
  /**
   * @brief             remove given object from
   *                        1.  list of existing oclObjects
   *                        2.  list of loaded oclObjects
   *                    delete buffer, if needed
   */
  void
  oclConnection ::
  removeDataObject      (const oclDataObject * const ocl_obj)
  {
    
    print_optional ("oclConnection :: removeDataObject (...)", m_verbose);
    
    // does ocl_obj has a buffer to be handled?
    if (ocl_obj -> getMemState ())
    {
        
      // decrease buffer reference count
      std::map <clBuffer *, int> :: iterator tmp_it =  m_current_buffers.find (ocl_obj -> getBuffer ());
      
      if (tmp_it == m_current_buffers.end ())
      {
      
        throw oclError ("Buffer not found!", "oclConnection :: removeDataObject");

      }
    
      print_optional (" => handle buffer!!! (ref_count: %d)", tmp_it -> second, m_verbose);
      
      // delete buffer object in case of no left references
      if (-- tmp_it -> second == 0)
      {
    
        print_optional (" *!* delete buffer *!*", m_verbose);
      
        // remove buffer from list
        m_current_buffers.erase (tmp_it);
        
        // delete buffer object
        delete tmp_it -> first;
    
      }
    
      // remove from list of loaded oclObjects
      m_loaded_ocl_objects.remove (ocl_obj -> getID ());
  
    }
    
    // remove given data object
    m_current_ocl_objects.erase (ocl_obj -> getID ());
    
  }
  
  
  /**
   * @brief             remove given object from
   *                        1.  list of loaded oclObjects
   *                    delete buffer, if needed
   */
  void
  oclConnection ::
  releaseBuffer        (const oclDataObject * const ocl_obj)
  {
    
    print_optional ("oclConnection :: releaseBuffer (...)", m_verbose);
    
    // does ocl_obj has a buffer to be handled?
    if (ocl_obj -> getMemState ())
    {
        
      // decrease buffer reference count
      std::map <clBuffer *, int> :: iterator tmp_it =  m_current_buffers.find (ocl_obj -> getBuffer ());
      
      if (tmp_it == m_current_buffers.end ())
      {
      
        throw oclError ("Buffer not found!", "oclConnection :: releaseBuffer");

      }
    
      print_optional (" => handle buffer!!! (ref_count: %d)", tmp_it -> second, m_verbose);
      
      // delete buffer object in case of no left references
      if (-- tmp_it -> second == 0)
      {
    
        print_optional (" *!* delete buffer *!*", m_verbose);
      
        // remove buffer from list
        m_current_buffers.erase (tmp_it);
        
        // delete buffer object
        delete tmp_it -> first;
    
      }
    
      // remove from list of loaded oclObjects
      m_loaded_ocl_objects.remove (ocl_obj -> getID ());
  
    }
    
  }
  
  
  // create function object for running an OpenCL kernel
  template <class T, class S>
  oclFunctionObject * const
  oclConnection ::
  makeFunctionObject    (const        string &               kernel_name,
                               oclDataObject * const * const        args, const      int &  num_args,
                         const    KernelType                 kernel_type, const SyncType   sync_type)
  {
    
    oclFunctionObject * kernel_obj;
   
    if (kernel_type == KERNEL)
    {
    
      if (sync_type == SYNC)
      {
        kernel_obj = new oclKernelObject (kernel_name, args, num_args);
      }
      else
      {
        kernel_obj = new oclAsyncKernelObject (kernel_name, args, num_args);
      }

    }
    
    return kernel_obj;
    
  }
  
  
  // create function object for running an ViennaCL algorithm
  template <class T, class S>
  oclFunctionObject * const
  oclConnection ::
  makeFunctionObject    (const   vclAlgoType &                 algo,
                               oclDataObject * const * const   args,
                         const           int &                 num_args,
                         const    KernelType                   kernel_type,
                         const      SyncType                   sync_type,
                         const           int                   num_scalars,
                                         int *                 scalars)
  {
  
    oclFunctionObject * algo_obj;
    
    if (kernel_type == VCL)
    {
    
      if (sync_type == SYNC)
      {
        algo_obj = new oclViennaClObject <T, S> (algo, args, num_args, num_scalars, scalars);
      }
      else
      {
  //      kernel_obj = new oclAsyncVCLObject (algo, args, num_args);
      }    
    
    }

    return algo_obj;
    
  }


  // create function object for running an ViennaCL algorithm
  template <class T, class S>
  oclFunctionObject * const
  oclConnection ::
  makeFunctionObject    (const oclAMDBlasType &                 algo,
                                oclDataObject * const * const   args,
                         const            int &                 num_args,
                         const     KernelType                   kernel_type,
                         const       SyncType                   sync_type,
                         const            int                   num_scalars,
                                          int *                 scalars)
  {
  
    oclFunctionObject * algo_obj;
    
    if (kernel_type == AMD)
    {
    
      if (sync_type == SYNC)
      {
        algo_obj = new oclAMDBlasObject <T, S> (algo, args, num_args, num_scalars, scalars);
      }
      else
      {
  //      kernel_obj = new oclAsyncVCLObject (algo, args, num_args);
      }    
    
    }

    return algo_obj;
    
  }
  
  
  /**
   * @brief             create buffer for given object id
   *
   * @param             size      ... in bytes
   */
  template <class T>
  void
  oclConnection ::
  createBuffer          (                T   * const cpu_arg,
                         const    ::size_t           size,
                         const oclObjectID           obj_id)
  {
    
    print_optional (" * oclConnection :: createBuffer (id: %d)", obj_id, m_verbose);
    
    // find corresponding data object
    std::map <oclObjectID, oclDataObject * const> :: iterator tmp_it = m_current_ocl_objects.find (obj_id);
  
    // check if data object exists
    if (tmp_it == m_current_ocl_objects.end ())
    {

      stringstream msg;
      msg << "oclDataObject (" << obj_id << ") does not exist!";

      throw oclError (msg.str (), "oclConnection :: createBuffer");

    }
      
    // retrieve data object
    oclDataObject * const p_tmp_obj = tmp_it -> second;
        
    // check if buffer exists
    if (p_tmp_obj -> getMemState ())
    {
    
      stringstream msg;
      msg << "oclDataObject (" << obj_id << ") already has a GPU buffer!";
      
      throw oclError (msg.str (), "oclConnection :: createBuffer");
 
    }

    // create buffer object
    p_tmp_obj -> setBuffer (new clBuffer (m_cont, CL_MEM_READ_WRITE, size, cpu_arg, &m_error));
    print_optional (" -!-> createBuffer: ", errorString (m_error), m_verbose);
      
    // add to buffer list
    m_current_buffers.insert (std::pair <clBuffer *, int> (p_tmp_obj -> getBuffer (), 1));
        
    // add data object to list of loaded objects
    m_loaded_ocl_objects.push_back (obj_id);
      
  }
  
  
  /**
   * @brief             load cpu data to GPU and return corresponding buffer
   *                      - TODO (just copied from setKernelArg (...)
   *
   * @param             size    ... in bytes
   */
  template <class T>
  void
  oclConnection ::
  loadToGPU             (       T   * const cpu_arg,
                         ::size_t           size,
                         clBuffer   * const buffer)
  {
    
    print_optional (" * oclConnection :: loadToGPU", m_verbose );
    
    // loop over command queues, write information to device (global memory)
    for (clCommandQueues::iterator it = m_comqs.begin(); it < m_comqs.end(); ++it)
    {
                   
      try {

        m_error = it -> enqueueWriteBuffer (*buffer, CL_FALSE, 0, size, cpu_arg, NULL, NULL);

      } catch (cl::Error cle) {

        throw oclError (cle.err (), "oclConnection :: loadToGPU");

      }
 
    }
     
  }



# endif // __OCL_CONNECTION_HPP__
