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
make_double_kernel        (std::string const & source, std::string const & fp_extension);








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
      
      template <class T>
      oclFunctionObject * const
      makeFunctionObject    (const string & kernel_name,
                             oclDataObject * const * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type);
                        
      template <class T>     
      oclFunctionObject * const
      makeFunctionObject    (const vclAlgoType & algo_name,
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
      template <class T>
      const clProgram       getProgram        () const { return ocl_precision_trait <T> :: getProgram (this); }
      
      
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
      template <class T>
      void
      activate              ( );


    
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
                               m_prog_f,     // program containing whole source code (single precision)
                               m_prog_d;     // program containing whole source code (double precision)
      clKernels              * mp_kernels,
                               m_kernels_f,  // kernels available in m_prog_f
                               m_kernels_d;  // kernels available in m_prog_d
      std::vector<clBuffers> * mp_buffers,
                               m_buffers_f,  // vectors of buffers for each kernels arguments (single prec)
                               m_buffers_d;  // vectors of buffers for each kernels arguments (double prec)
      cl_int                   m_error;      // error variable
      VerbosityLevel           m_verbose;    // verbosity level
      clKernel               * mp_actKernel;
      int                      num_kernel;

      
      template <class T>
      struct ocl_precision_trait
      { };


      template <class T> //, class trait = read_kernel_source_trait <T> >
      static const char * ReadSource   (const char *, int * size);
             int          BuildProgram ();
             const char * errorString (cl_int e);

  
      /******************
       ** constructors **
       ******************/
      oclConnection         (const char *,
                             cl_device_type device_type = CL_DEVICE_TYPE_GPU,
                             VerbosityLevel verbose = VERB_NONE);
      
      oclConnection         (oclConnection &) { /*TODO*/};
      

  }; // class oclConnection
  
  
  
  
  
      template <>
      struct oclConnection :: ocl_precision_trait <float>
      {

        public:
  

          typedef float elem_type;

          static inline
          const char *
          modify_source (void * buf, int * size)
          {
            return (const char *) buf;
          }
          
          static inline
          clProgram *
          getProgram    (oclConnection * const con)
          {
            return & (con -> m_prog_f);
          }
          
          static inline
          std::vector <clBuffers> *
          getBuffers    (oclConnection * const con)
          {
            return & (con -> m_buffers_f);
          }
          
          static inline
          clKernels *
          getKernels    (oclConnection * const con)
          {
            return & (con -> m_kernels_f);
          }

      };

      template <>
      struct oclConnection :: ocl_precision_trait <double>
      {

        public:
  

          typedef double elem_type;

          static inline
          const char *
          modify_source (void * buf, int * size)
          {
            std::string * tmp_str = new string (make_double_kernel (std::string ((const char *) buf), std::string ("cl_khr_fp64")));
            *size = tmp_str -> size () * sizeof (char);
            return tmp_str->c_str ();
          }
          
          static inline
          clProgram *
          getProgram    (oclConnection * const con)
          {
            return & (con -> m_prog_d);
          }
          
          static inline
          std::vector <clBuffers> *
          getBuffers    (oclConnection * const con)
          {
            return & (con -> m_buffers_d);
          }
          
          static inline
          clKernels *
          getKernels    (oclConnection * const con)
          {
            return & (con -> m_kernels_d);
          }
          
      };

      template <>
      struct oclConnection :: ocl_precision_trait <bool>
      {

        public:
  

          typedef bool elem_type;

          static inline
          const char *
          modify_source (void * buf, int * size)
          {
            return (const char *) buf;
          }
          
          static inline
          clProgram *
          getProgram    (oclConnection * const con)
          {
            return & (con -> m_prog_f);
          }
          
          static inline
          std::vector <clBuffers> *
          getBuffers    (oclConnection * const con)
          {
            return & (con -> m_buffers_f);
          }
          
          static inline
          clKernels *
          getKernels    (oclConnection * const con)
          {
            return & (con -> m_kernels_f);
          }

      };
  
  
  
  

  /*******************
   ** includes (II) **
   *******************/
   
  // ocl
  # include "oclFunctionObject.hpp"
  # include "oclDataObject.hpp"
  # include "oclKernelObject.hpp"
  # include "oclViennaClObject.hpp"
  # include "oclAsyncKernelObject.hpp"
  
  

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



  // set central vars for particular precision mode for type T
  template <class T>
  void
  oclConnection ::
  activate              ( )
  {
    mp_prog    = ocl_precision_trait <T> :: getProgram (this);
    mp_kernels = ocl_precision_trait <T> :: getKernels (this);
    mp_buffers = ocl_precision_trait <T> :: getBuffers (this);
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
      mp_inst = new oclConnection ( "./src/matrix/ocl/test.cl", CL_DEVICE_TYPE_GPU, VERB_MIDDLE );
      
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
  template <class T>
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
  template <class T>
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
        algo_obj = new oclViennaClObject <T> (algo, args, num_args, num_scalars, scalars);
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
