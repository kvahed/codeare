/**************
 ** includes **
 **************/

// ocl
# include "oclConnection.hpp"



/**************************
 ** function definitions **
 **************************/


/** 
 * @brief                  Create a double precision kernel out of a single precision kernel.
 *
 * @param  source          The source string.
 * @param  fp_extension    An info string that specifies the OpenCL double precision extension.
 * @param  vector_type     Type of vector elements in kernel file.
 * @param  scalar_type     Type of scalar elements in kernel file.
 * @param  n               multiple of built-in OpenCl vector types in kernel file.
 *
 * @return                 The double precision kernel.
 */
inline
std::string
modify_kernel              ( std::string const &       source,
                             std::string const & fp_extension, 
                             std::string const &  A_type,
                             std::string const &  B_type,
                                     int const &            n )
{

  std::stringstream ss;          
  ss <<              "# pragma OPENCL EXTENSION " << fp_extension << "\n\n";
  ss <<              "# define  A_type_n " << A_type;
  if (n > 0) ss << n;
  ss << std::endl << "# define  A_type "   << A_type <<      std::endl;
  ss <<              "# define  B_type_n " << B_type;
  if (n > 0) ss << n;
  ss << std::endl << "# define  B_type "   << B_type <<      std::endl;
  ss <<              "# define  vec_len "  << (n > 0 ? n : 1) << std::endl;
  ss << source << std::endl;
  
  return ss.str ();

}



template <class T, class S>
const char*
oclConnection::
ReadSource            (std::string   fname,
                               int * size  )
{

  // open file
  FILE *f = fopen(fname.c_str (), "r");

  // create buffer for source code
  void *buf;

  if (!f) {
    fprintf(stderr, "Unable to open %s for reading\n", fname.c_str ());
    return NULL;
  }

  // explore file information
  fseek (f, 0, SEEK_END);
  *size = ftell (f);
  fseek (f, 0, SEEK_SET);

  // read and close file
  buf = malloc (*size+1);
  *size = fread (buf, 1, *size, f);
  fclose (f);
  ((char*)buf)[*size] = '\0';

  /* source code !and! size may change with different precision */
  return ocl_precision_trait <T, S> :: modify_source (buf, size); 

}



template <class T, class S>
oclError &
oclConnection ::
init_program_kernels    (oclConnection * const con)
{

  print_optional (" ** Initializing oclConnection for types ", ocl_precision_trait <T, S> :: getTypeString (), con -> m_verbose);

  // get sources
  int size;
  std::vector <std::pair <const char *, ::size_t> > sources;
  const std::vector <std::string> filenames = ocl_precision_trait <T, S> :: getFilenames ();
  for (std::vector <std::string> ::const_iterator it = filenames.begin (); it != filenames.end (); ++it)
  {
    const char * tmp_src = ReadSource <T, S> (*it, &size);
    sources.push_back (std::pair <const char*, ::size_t> (tmp_src, size));
  }
  
  clProgram prog;
  clKernels kernels;
  
  try
  {

    // create program
    prog = clProgram (con -> m_cont, sources, & con -> m_error);
    *(ocl_precision_trait <T, S> :: getProgram (con)) = prog;
  
    // build program
    con -> m_error = prog.build (con -> m_devs);

    // create kernels
    con -> m_error = prog.createKernels (&kernels);
    *(ocl_precision_trait <T, S> :: getKernels (con)) = kernels;

  }
  catch (cl::Error cle)
  {
    cout << " Type: " << ocl_precision_trait <T, S> :: getTypeString () << std::endl;
    cout << "Error while building program: " << cle.what ()                                            << endl;
    cout << "Build Status: "                 << prog.getBuildInfo<CL_PROGRAM_BUILD_STATUS>  (con -> m_devs [0]) << endl;
    cout << "Build Options:\t"               << prog.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS> (con -> m_devs [0]) << endl;
    cout << "Build Log:\t "                  << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>     (con -> m_devs [0]) << endl;
    std::cout << " Error flag: " << cle.err () << " (" << con -> errorString (cle.err ()) << ")" << std::endl;
    throw -1;
  }
  
}




// constructor
oclConnection::
oclConnection ( const char      * filename_A_type,
                const char      * filename_AB_type,
                const char      * filename_A_cx_type,
                const char      * filename_AB_cx_type,
                const char      * filename_A_mx_type,
                const char      * filename_AB_mx_type,
                cl_device_type    device_type,
                VerbosityLevel    verbose )
              : m_current_ocl_objects (),
                m_current_buffers (),
                m_loaded_ocl_objects ()
{

  m_verbose = verbose;

  // platform
  clPlatforms tmp_platforms;
  m_error = clPlatform::get (& tmp_platforms);
  print_optional (" ** # of platforms: %d", tmp_platforms.size(), VERB_LOW);
  if (tmp_platforms.size() > 0)
    m_plat = tmp_platforms [0];    // choose first available device
  else
    throw oclError ("No platform available", "oclConnection :: CTOR");

  // devices
  m_error = m_plat.getDevices (device_type, &m_devs);
  print_optional (" ** # of devices on platform: %d", m_devs.size(), VERB_LOW);
  std::string vendor;
  print_optional (" ** device type (0): ", (m_devs [0].getInfo (CL_DEVICE_VENDOR, &vendor), vendor.c_str ()), VERB_LOW);
  if (m_devs.size() == 0)
    throw oclError ("No devices available on this platform", "oclConnection :: CTOR");

  // context /** ViennaCL **/ /* TODO */
  m_cont = clContext ( viennacl::ocl::current_context () . handle () . get ()); //clContext (m_devs);       // same context for all devices
  
  // command queues
  for (clDevices::iterator it = m_devs.begin(); it < m_devs.end(); ++it)  // iterate over all devices and create a command queue for each
  {
    m_comqs.push_back (clCommandQueue (m_cont, *it));
  }

  init_program_kernels < float,  float> (this);
  init_program_kernels <double, double> (this);
  init_program_kernels <  cxfl,  float> (this);
  init_program_kernels <  cxfl,   cxfl> (this);
  init_program_kernels <  cxdb, double> (this);
  init_program_kernels <  cxdb,   cxdb> (this);

  // set current kernel!
  num_kernel = 0;
  mp_actKernel = & (m_kernels_ff [num_kernel]); // defined default value !!!

  print_optional (" ** oclConnection constructed!", m_verbose);
    
  /**
   * setup ViennaCL
   */
  print_optional (" ** setup ViennaCl!", m_verbose);
  viennacl :: vector <float> tmp (10);
  viennacl :: vector <double> tmp2 (10);
  tmp = tmp + tmp;

}




// activate kernel, stays active until another kernel is activated (default: 0)
//template <class T>
int
oclConnection::
activateKernel            (const std::string kernelname)
{

  std::string act_kernelname;
  
  // check if kernel already activated
//  m_error = mp_actKernel -> getInfo <std::string> (CL_KERNEL_FUNCTION_NAME, &act_kernelname);

  /************************/
  // ---> removed check by name, if kernel is already activated,
  //      since the same kernels for different precisions have the same name
  /************************/
  
//  if (kernelname.compare (act_kernelname) == 0)
//    return 0;

  // search for matching kernel  
  for (int i = 0; i < mp_kernels -> size(); i++)
  {
    m_error = (* mp_kernels) [i] . getInfo <std::string> (CL_KERNEL_FUNCTION_NAME, &act_kernelname);
    if (kernelname.compare (act_kernelname) == 0)
    {
      mp_actKernel = & ((* mp_kernels) [i]);
      return 0;
    }
  }

  /* throw error message */
  stringstream msg;
  msg << "No kernel found (" << kernelname << ")";
  throw oclError (msg.str (), "oclConnection :: activateKernel");

  return -1;

}




// run active kernel with given dimensions for work groups and work items
int
oclConnection::
runKernel             (const cl::NDRange  & global_dims,
                       const cl::NDRange  & local_dims  )
{

  // execute activated kernel on all available devices 
  for (clCommandQueues::iterator it = m_comqs.begin(); it < m_comqs.end(); ++it)
  {
    try {
    
      m_error = (*it).enqueueNDRangeKernel (*mp_actKernel, cl::NullRange, global_dims, local_dims);

    } catch (cl::Error cle) {

      throw oclError (cle.err (), "oclConnection :: runKernel");

    }
  }
  
}



const char*
oclConnection::
errorString             (cl_int e)
{

  static const char* errorString[] = {
    "CL_SUCCESS",
    "CL_DEVICE_NOT_FOUND",
    "CL_DEVICE_NOT_AVAILABLE",
    "CL_COMPILER_NOT_AVAILABLE",
    "CL_MEM_OBJECT_ALLOCATION_FAILURE",
    "CL_OUT_OF_RESOURCES",
    "CL_OUT_OF_HOST_MEMORY",
    "CL_PROFILING_INFO_NOT_AVAILABLE",
    "CL_MEM_COPY_OVERLAP",
    "CL_IMAGE_FORMAT_MISMATCH",
    "CL_IMAGE_FORMAT_NOT_SUPPORTED",
    "CL_BUILD_PROGRAM_FAILURE",
    "CL_MAP_FAILURE",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "CL_INVALID_VALUE",
    "CL_INVALID_DEVICE_TYPE",
    "CL_INVALID_PLATFORM",
    "CL_INVALID_DEVICE",
    "CL_INVALID_CONTEXT",
    "CL_INVALID_QUEUE_PROPERTIES",
    "CL_INVALID_COMMAND_QUEUE",
    "CL_INVALID_HOST_PTR",
    "CL_INVALID_MEM_OBJECT",
    "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
    "CL_INVALID_IMAGE_SIZE",
    "CL_INVALID_SAMPLER",
    "CL_INVALID_BINARY",
    "CL_INVALID_BUILD_OPTIONS",
    "CL_INVALID_PROGRAM",
    "CL_INVALID_PROGRAM_EXECUTABLE",
    "CL_INVALID_KERNEL_NAME",
    "CL_INVALID_KERNEL_DEFINITION",
    "CL_INVALID_KERNEL",
    "CL_INVALID_ARG_INDEX",
    "CL_INVALID_ARG_VALUE",
    "CL_INVALID_ARG_SIZE",
    "CL_INVALID_KERNEL_ARGS",
    "CL_INVALID_WORK_DIMENSION",
    "CL_INVALID_WORK_GROUP_SIZE",
    "CL_INVALID_WORK_ITEM_SIZE",
    "CL_INVALID_GLOBAL_OFFSET",
    "CL_INVALID_EVENT_WAIT_LIST",
    "CL_INVALID_EVENT",
    "CL_INVALID_OPERATION",
    "CL_INVALID_GL_OBJECT",
    "CL_INVALID_BUFFER_SIZE",
    "CL_INVALID_MIP_LEVEL",
    "CL_INVALID_GLOBAL_WORK_SIZE",
  };

  const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

  const int index = -e;

  return (index >= 0 && index < errorCount) ? errorString[index] : "";

}







