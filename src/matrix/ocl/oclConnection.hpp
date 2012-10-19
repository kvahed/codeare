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
      static
      oclConnection *
      Instance              ();

      void
      Initialise            ();
    
      void
      Finalise              ();
      
      oclFunctionObject * const
      makeFunctionObject    (const string & kernel_name,
                             oclDataObject * const * const args, const int & num_args,
                             const KernelType kernel_type, const SyncType sync_type);
                              
      void
      addDataObject         (oclDataObject * const ocl_obj, const oclObjectID & obj_id);
                          
      void
      run                   (oclFunctionObject * const func_obj) const;

      // getters for OpenCL objects
      const clPlatform      getPlatform       () const { return m_plat; }
      const clDevices       getDevices        () const { return m_devs; }
      const clContext       getContext        () const { return m_cont; }
      const clCommandQueues getCommandQueues  () const { return m_comqs; }
      const clProgram       getProgram        () const { return m_prog; }
      
      
      /****************************************
       ** functions for OpenCL functionality **
       ****************************************/
  
      // activate kernel (stays activated until another kernel is activated
      int activateKernel      (const std::string kernelname);
    
      // run kernel with given dimensions
      int runKernel           (const cl::NDRange & global_dims,
                               const cl::NDRange & local_dims);
    
    
      // set argument for kernel by pointer, optionally return created buffer
      template <typename T>
      cl_int
      setKernelArg            (int           num,
                               T           * in_argument,
                               ::size_t      size,
                               cl_mem_flags  mem_flag,
                               clBuffer    * out = NULL)
      {
    
        // create new buffer
        clBuffer * p_tmp_buffer = new clBuffer(m_cont, mem_flag, size, in_argument, &m_error);

        // if buffer is used for input, copy data to device 
        if (   mem_flag == CL_MEM_READ_ONLY
            || mem_flag == CL_MEM_READ_WRITE)
        {
      
          // loop over command queues, write information to device (global memory)
          for (clCommandQueues::iterator it = m_comqs.begin(); it < m_comqs.end(); ++it)
          {
        
            m_error = it -> enqueueWriteBuffer (*p_tmp_buffer, CL_TRUE, 0, size, in_argument, NULL, NULL);

          }
        
        }
      
        // register buffer as argument      
        m_error = mp_actKernel -> setArg (num, *p_tmp_buffer);
        m_buffers [num_kernel].push_back (p_tmp_buffer);

        // if buffer is used for output, return buffer object
        if (   mem_flag == CL_MEM_READ_WRITE
            || mem_flag == CL_MEM_WRITE_ONLY)
        {

          out = p_tmp_buffer;

        }
        else // delete buffer, if not needed
        {

          out = NULL;
  
        }
  
      }
  
  
      // TODO: distinction between read/write, read_only, write_only !!!
      // set argument for kernel by buffer
      template <typename T>
      cl_int
      setKernelArg            (int            num,
                               clBuffer     * buf,
                               T            * in_argument,
                               ::size_t       size)
      {
      
        cout << " ** !! not to use !! ** " << endl;

        for (clCommandQueues::iterator it = m_comqs.begin(); it < m_comqs.end(); ++it)
        {

          m_error = it -> enqueueWriteBuffer (*buf, CL_TRUE, 0, size, in_argument, NULL, NULL);

        }

        m_error = mp_actKernel -> setArg (num, *buf);

      }
    
    
      // get argument from activated kernel
      template <typename T>
      cl_int
      getKernelArg            (int        num,
                               T        * out_argument,
                               ::size_t   size)
      {

        try {

          m_error = m_comqs [0] . enqueueReadBuffer (*(m_buffers [num_kernel] [num]), CL_TRUE, 0, size, out_argument, NULL, NULL);

        } catch (cl::Error cle) {

          cout << "error: " << errorString (cle.err()) << endl;

        }

        return m_error;

      }
    
    
    private:
    
      /**********************
       ** member variables **
       **********************/
      static oclConnection                    * mp_inst;
      std::map <oclObjectID, oclDataObject *>   m_current_ocl_objects;
      std::list <oclObjectID>                   m_loaded_ocl_objects;
      
      /********************
       ** OpenCL members **
       ********************/
      clPlatform              m_plat;       // available platforms
      clDevices               m_devs;       // vector of available devices on m_plat
      clContext               m_cont;       // context used associated with m_devs
      clCommandQueues         m_comqs;      // each command queue is associated with corresponding device from m_devs
      clProgram               m_prog;       // program containing whole source code
      clKernels               m_kernels;    // kernels available in m_prog
      std::vector<clBuffers>  m_buffers; // vectors of buffers for each kernels arguments
      cl_int                  m_error;      // error variable
      bool                    m_verbose;    // verbosity level
      clKernel              * mp_actKernel;
      int                     num_kernel;

      static const char * ReadSource   (const char *, int * size);
             int          BuildProgram ();
             const char * errorString (cl_int e);

  
      /******************
       ** constructors **
       ******************/
      oclConnection         (const char *,
                             cl_device_type device_type = CL_DEVICE_TYPE_GPU,
                             bool verbose = false);
      
      oclConnection         (oclConnection &) { /*TODO*/};
        

  }; // class oclConnection
  
  

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
  

  
  oclConnection *
  oclConnection ::
  Instance ()
  {
  
    if (mp_inst == NULL)
      mp_inst = new oclConnection ( "./src/matrix/ocl/test.cl", CL_DEVICE_TYPE_GPU, true );
      
    return mp_inst;
    
  }
  
  
  void
  oclConnection ::
  addDataObject         (oclDataObject * const ocl_obj, const oclObjectID & obj_id)
  {
    std::cout << "addDataObject" << std::endl;
    m_current_ocl_objects.insert (std::pair <oclObjectID, oclDataObject *> (obj_id, ocl_obj));
  }
  
  
  // create function object for running an OpenCL kernel
  oclFunctionObject * const
  oclConnection ::
  makeFunctionObject    (const string & algo_name,
                         oclDataObject * const * const args, const int & num_args,
                         const KernelType kernel_type, const SyncType sync_type)
  {
    
    oclFunctionObject * kernel_obj;
   
    if (kernel_type == KERNEL)
    {
    
      if (sync_type == SYNC)
      {
        kernel_obj = new oclKernelObject (algo_name, args, num_args);
      }
      else
      {
  //      kernel_obj = new oclAsyncKernelObject (kernel_name, args, num_args);
      }

    }
    else
    {
    
      if (sync_type == SYNC)
      {
        kernel_obj = new oclViennaClObject (algo_name, args, num_args);
      }
      else
      {
  //      kernel_obj = new oclAsyncVCLObject (algo, args, num_args);
      }    
    
    }

    return kernel_obj;
    
  }



# endif // __OCL_CONNECTION_HPP__
