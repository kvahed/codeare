/**************
 ** includes **
 **************/


# ifdef __USE_VIENNA_CL__
  // ViennaCL
  # include <viennacl/ocl/backend.hpp>
# endif

#include "oclConnection.hpp"


/**************************
 ** function definitions **
 **************************/


/** 
 * @brief                  Modify kernel source code by adding predefined makro definitions.
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


/**
 * @brief                 Add arbitrary makro definition at the front of the
 *                        given kernel source code.
 * 
 * @param  source         Kernel source code.
 * @param  vec_makros     Vector of makro definitions (without keyword).
 * 
 * @return                Modified kernel source code.
 */
std::string
modify_kernel             ( const std::string               & source,
                            const std::vector <std::string> & vec_makros )
{
  
  std::stringstream ss;
  for (std::vector <std::string> :: const_iterator it = vec_makros.begin ();
          it != vec_makros.end (); it++)
    ss << "# define " << *it << std::endl;
  
  ss << source << std::endl;
  
  return ss.str ();
  
}



template <class T, class S>
const char*
oclConnection::
ReadSource            (std::string   fname,
                               int *  size)
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
  const char * tmp_src = ocl_precision_trait <T, S> :: modify_source (buf, size); 
  
  free (buf);
  
  return tmp_src;

}


template <class T, class S>
oclError &
oclConnection ::
rebuildWithSource       (const std::string & filename)
{
  init_program_kernels <T, S> (this, std::vector <std::string> (1, filename));
}


template <class T, class S>
oclError &
oclConnection ::
rebuildWithSource       (const std::string & filename,
                         const std::vector <std::string> & makros)
{
  init_program_kernels <T, S> (this, std::vector <std::string> (1, filename), makros);
}


template <class T, class S>
oclError &
oclConnection ::
rebuildWithSources      (const std::vector <std::string> & filenames)
{
  init_program_kernels <T, S> (this, filenames);
}


template <class T, class S>
oclError &
oclConnection ::
rebuildWithSources      (const std::vector <std::string> & filenames,
                         const std::vector <std::string> & makros)
{
  init_program_kernels <T, S> (this, filenames, makros);
}


template <class T, class S>
oclError &
oclConnection ::
init_program_kernels    (oclConnection * const con, const std::vector <std::string> & add_filenames, const std::vector <std::string> & add_makros)
{

  print_optional (" ** Initializing oclConnection for types ", ocl_precision_trait <T, S> :: getTypeString (), con -> m_verbose);

  // get sources
  int size;
  std::vector <std::pair <const char *, ::size_t> > sources;

  // read standard sources from precision trait
  const std::vector <std::string> filenames = ocl_precision_trait <T, S> :: getFilenames ();
  for (std::vector <std::string> ::const_iterator it = filenames.begin (); it != filenames.end (); ++it)
  {
    const char * tmp_src = ReadSource <T, S> (*it, &size);
    sources.push_back (std::pair <const char*, ::size_t> (tmp_src, size));
  }
  
  // read additional sources
  if (add_filenames.size ())
  {
    for (std::vector <std::string> :: const_iterator it = add_filenames.begin (); it != add_filenames.end (); ++it)
    {
      const char * tmp_src = ReadSource <T, S> (*it, &size);
      if (add_makros.size ())
      {
        std::string * p_tmp_src = new std::string (modify_kernel(std::string (tmp_src), add_makros));
        size = p_tmp_src -> size () * sizeof (char);
        tmp_src = p_tmp_src -> c_str ();
      }
      sources.push_back (std::pair <const char *, ::size_t> (tmp_src, size));
    }
  }
  
  clProgram prog;
  clKernels kernels;
    
  try
  {

    // create program
    prog = clProgram (con -> m_cont, sources, & con -> m_error);

    // build program
	  con -> m_error = prog.build (con -> m_devs, "-cl-std=CL1.0");//1 -cl-mad-enable -cl-fast-relaxed-math -cl-nv-maxrregcount=32 -cl-nv-verbose");
    
//    std::cout << "Build Log:\t "                  << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>     (con -> m_devs [0]) << std::endl;

    // create kernels
    con -> m_error = prog.createKernels (&kernels);

    // assign created objects to type specific storage
    *(ocl_precision_trait <T, S> :: getProgram (con)) = prog;
    *(ocl_precision_trait <T, S> :: getKernels (con)) = kernels;

  }
  catch (cl::Error & cle)
  {
    std::cout << " Type: " << ocl_precision_trait <T, S> :: getTypeString () << std::endl;
    std::cout << "Error while building program: " << cle.what ()                                                     << std::endl;
    std::cout << "Build Status: "                 << prog.getBuildInfo<CL_PROGRAM_BUILD_STATUS>  (con -> m_devs [0]) << std::endl;
    std::cout << "Build Options:\t"               << prog.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS> (con -> m_devs [0]) << std::endl;
    std::cout << "Build Log:\t "                  << prog.getBuildInfo<CL_PROGRAM_BUILD_LOG>     (con -> m_devs [0]) << std::endl;
    std::cout << " Error flag: " << cle.err () << " (" << con -> errorString (cle.err ()) << ")" << std::endl;
    throw -1;
  }
    
}




// constructor
oclConnection::
oclConnection ( cl_device_type    device_type,
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
  try {
	  m_error = m_plat.getDevices (device_type, &m_devs);
  } catch (cl::Error & cle)
  {
	  std::cerr << oclError (cle.err(), " Error while requesting device ids! ") << std::endl;
  	  std::cerr << cle.what () << std::endl;
  }
  print_optional (" ** # of devices on platform: %d", m_devs.size(), VERB_LOW);
  std::string info_string;
  print_optional (" ** device type [0]: ", (m_devs [0].getInfo (CL_DEVICE_VENDOR, &info_string), info_string.c_str ()), VERB_LOW);
  print_optional (" ** device extensions [0]: ", (m_devs [0].getInfo (CL_DEVICE_EXTENSIONS, &info_string), info_string.c_str ()), VERB_MIDDLE);
  print_optional (" ** device max WG size [0]: %d", (m_devs [0].getInfo (CL_DEVICE_MAX_WORK_GROUP_SIZE, &m_max_wg_size), m_max_wg_size), VERB_NONE);
  print_optional (" ** device max WI dim [0]: %d", (m_devs [0].getInfo (CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &m_max_wi_dim), m_max_wi_dim), VERB_NONE);
  m_devs [0].getInfo (CL_DEVICE_MAX_WORK_ITEM_SIZES, &m_max_wi_sizes);
  print_optional (" ** device max WI sizes [0]: %zu, %zu, %zu", m_max_wi_sizes[0], m_max_wi_sizes[1], m_max_wi_sizes[2], VERB_NONE);
  if (m_devs.size() == 0)
    throw oclError ("No devices available on this platform", "oclConnection :: CTOR");
  
  // don't use multiple devices (for now)
  if (m_devs.size () > 1)
  {
    clDevices tmp_devs;
    tmp_devs.push_back (m_devs [0]);
    m_devs = tmp_devs;
  }
  
# ifdef __USE_VIENNA_CL__
  // context /** ViennaCL **/ /* TODO */
  m_cont = clContext ( viennacl::ocl::current_context () . handle () . get ()); //clContext (m_devs);       // same context for all devices
# else
  m_cont = clContext (m_devs); // same context for all devices
# endif
  
  // command queues
  for (clDevices::iterator it = m_devs.begin(); it < m_devs.end(); ++it)  // iterate over all devices and create a command queue for each
  {
    m_comqs.push_back (clCommandQueue (m_cont, *it, CL_QUEUE_PROFILING_ENABLE));
  }

  init_program_kernels < float,  float> (this);
  init_program_kernels <double, double> (this);
//  init_program_kernels <  cxfl,  float> ();
  init_program_kernels <  cxfl,   cxfl> (this);
//  init_program_kernels <  cxdb, double> ();
  init_program_kernels <  cxdb,   cxdb> (this);

  // set current kernel!
  num_kernel = 0;
  mp_actKernel = & (m_kernels_ff [num_kernel]); // defined default value !!!

  print_optional (" ** oclConnection constructed!", m_verbose);
    
# ifdef __USE_VIENNA_CL__
  /**
   * setup ViennaCL
   */
  print_optional (" ** setup ViennaCl!", m_verbose);
  viennacl :: vector <float> tmp (10);
  viennacl :: vector <double> tmp2 (10);
  tmp = tmp + tmp;
# endif
  
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
  std::stringstream msg;
  msg << "No kernel found (" << kernelname << ")";
  throw oclError (msg.str (), "oclConnection :: activateKernel");

  return -1;

}


void
oclConnection::
waitForEvent          (const cl::Event & event)
const
{
  event.wait ();
}


const ProfilingInformation
oclConnection::
getProfilingInformation (const cl::Event & event)
const
{
  
  cl_ulong start, end;
  
  event.wait ();
  
  // read timing information
  clGetEventProfilingInfo(event(), CL_PROFILING_COMMAND_START, sizeof (cl_ulong), &start, NULL);
  clGetEventProfilingInfo(event(), CL_PROFILING_COMMAND_END, sizeof (cl_ulong), &end, NULL);
  
  return {start*1.e-9, end*1.e-9, 0};
  
}


// run active kernel with given dimensions for work groups and work items
const cl::Event
oclConnection::
runKernel             (const cl::NDRange  & global_dims,
                       const cl::NDRange  & local_dims  )
{
  
  // execute activated kernel on all available devices 
  for (clCommandQueues::iterator it = m_comqs.begin(); it < m_comqs.end(); ++it)
  {
    try {

      cl::Event event;
      
      m_error = (*it).enqueueNDRangeKernel (*mp_actKernel, cl::NullRange, global_dims, local_dims, NULL, &event);
      
      m_error = (*it).enqueueBarrier ();
      
      return event;
      
    } catch (cl::Error & cle) {

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







