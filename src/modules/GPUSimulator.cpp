#include "GPUSimulator.hpp"
using namespace RRStrategy;


GPUSimulator::GPUSimulator () {

	printf ("  Intialising GPU device & context ...\n");

	// Platform -----------------------

    std::vector<cl::Platform> pfs;
	m_error = cl::Platform::get(&pfs);
	
    printf("    number of platforms: %d\n", pfs.size());
	// --------------------------------

	// Use first available device -----

	m_dev = 0;

	cl_context_properties properties[] = 
        {CL_CONTEXT_PLATFORM, (cl_context_properties)(pfs[0])(), 0};
	// --------------------------------

	// Context and device -------------

    m_ctxt  = cl::Context(CL_DEVICE_TYPE_GPU, properties);
    m_devs  = m_ctxt.getInfo<CL_CONTEXT_DEVICES>();

    printf("    number of devices %d\n", m_devs.size());
	// --------------------------------

	// Command queue used for OpenCL commands

	try {
        m_cmdq = cl::CommandQueue(m_ctxt, m_devs[m_dev], 0, &m_error);
    } catch (cl::Error cle) {
        printf("  ERROR: %s(%d)\n", cle.what(), cle.err());
    }
	// --------------------------------

}
 

void 
GPUSimulator::ReadAndBuild (std::string ksrc) {
	
	int ksize = ksrc.size();

	printf("kernel size: %d\n", ksize);
	
	try {
		cl::Program::Sources cps (1, std::make_pair(ksrc.c_str(), ksize));
		m_prg = cl::Program(m_ctxt, cps);
	} catch (cl::Error cle) {
		printf("ERROR: %s(%s)\n", cle.what(), ErrorString (cle.err()));
	}
	
	printf("build program\n");

	try {
		m_error = m_prg.build(m_devs);
	} catch (cl::Error er) {
		printf("m_prg.build: %s\n", ErrorString (er.err()));
	}

	printf("done building program\n");

	std::cout << "Build Status: "   << m_prg.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(m_devs[0])  << std::endl;
	std::cout << "Build Options:\t" << m_prg.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(m_devs[0]) << std::endl;
	std::cout << "Build Log:\t "    << m_prg.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_devs[0])     << std::endl;
	
}



void
GPUSimulator::Simulate     (const Matrix<cplx>&   txm, const Matrix<cplx>&   rxm, 
							const Matrix<cplx>&    rf, const Matrix<double>&  gr, 
							const Matrix<double>&   r, const Matrix<double>&  m0, 
							const Matrix<double>& b0m, const double&          dt, 
							const bool&           exc, const bool&             v, 
							const size_t&          np, 
						          Matrix<cplx>&   res, Matrix<double>&         m) {

	ticks            tic  = getticks();

	size_t           nr   =   r.Dim(1);
	size_t           ntxc = txm.Dim(1);
	size_t           nrxc = rxm.Dim(1);
	size_t           nt   =  gr.Dim(1);
	
	printf ("  Simulaing: %s on %04i isochromats ... ", (exc) ? "         excitation" : " signal acquisition", (int)nr); fflush(stdout);
	
	// Anything to do? ----------------
	
	if (gr.Size() < 1 || r.Size() < 1) {
		std::cout << "  Bailing out: %i Gradient step for %i isochromats? I don't think so!\n" << std::endl;
		return;
	}
	// --------------------------------
	
	Matrix<cplx> mres;
	
	if (!exc)
		mres = Matrix<cplx> (res.Dim(0), res.Dim(1), np);
	
	printf (" done. (%.4f s)\n", elapsed(getticks(), tic) / Toolbox::Instance()->ClockRate());
	
}


char* 
GPUSimulator::ReadSource (const char* fname, int* size) {
	
	FILE *f = fopen(fname, "r");
	
	void *buf;
	
	if (!f) {
		fprintf(stderr, "Unable to open %s for reading\n", fname);
		return NULL;
	}
	
	fseek (f, 0, SEEK_END);
	*size = ftell (f);
	fseek (f, 0, SEEK_SET);
	
	buf = malloc (*size+1);
	*size = fread (buf, 1, *size, f);
	fclose (f);
	((char*)buf)[*size] = '\0';
	
	return (char*)buf;
	
}
		

const char*	
GPUSimulator::ErrorString (cl_int e) {
	
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
	
	const int index      = -m_error;
	
	return (index >= 0 && index < errorCount) ? errorString[index] : "";
	
}

