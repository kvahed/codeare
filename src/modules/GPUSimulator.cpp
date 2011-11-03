#include "GPUSimulator.hpp"

using namespace RRStrategy;



std::pair<const char*, unsigned> ReadCLFromFile (std::string fname) {

	std::string str, src;
	std::ifstream in;
	
	if(!in)
		std::cerr << "File not found??" << std::endl;
	
	in.open (fname.c_str());
	std::getline (in, str);
	
	while (in) {
		
		src += str + "\n";
		std::getline (in, str);
		
	}
	
	in.close();
	
	return std::make_pair(src.c_str(), src.length());

}


GPUSimulator::GPUSimulator (SimulationBundle* sb) {

	m_sb = sb;

	m_gdt = GAMMARAD * m_sb->dt;
	m_nt  = m_sb->agr->Dim(1);      // Time points
	m_nc  = m_sb->tb1->Dim(1);     // # channels
	m_ne  = m_sb->sr->Dim(1); 
	m_na  = m_sb->tr->Dim(1);
	m_nl  = 16;

	printf ("  Intialising GPU device & context ...\n");

	// Platform -----------------------

    std::vector<cl::Platform> pfs;
	m_error = cl::Platform::get(&pfs);
	
    printf("    number of platforms: %d\n", (int) pfs.size());
	// --------------------------------

	// Use first available device -----

	cl_context_properties properties[] = 
        {CL_CONTEXT_PLATFORM, (cl_context_properties)(pfs[0])(), 0};
	// --------------------------------

	// Context and device -------------

    m_ctxt  = cl::Context(CL_DEVICE_TYPE_GPU, properties);
    m_devs  = m_ctxt.getInfo<CL_CONTEXT_DEVICES>();

    printf("    number of devices %d\n", (int) m_devs.size());
	// --------------------------------

	// Command queue used for OpenCL commands

	try {
		m_cmdq = cl::CommandQueue(m_ctxt, m_devs[0], 0, &m_error);
    } catch (cl::Error cle) {
        printf("  ERROR: %s(%d)\n", cle.what(), cle.err());
    }
	// --------------------------------

}
 

bool
GPUSimulator::BuildProgram (std::string ksrc) {
	
	std::pair<const char*, unsigned> kpair = ReadCLFromFile (ksrc);

	printf("    OpenCL kernel size: %d\n    Assembling program ... ", (int) kpair.second); fflush (stdout);
	
	try {
		cl::Program::Sources cps (1, kpair);
		m_prg = cl::Program(m_ctxt, cps);
	} catch (cl::Error cle) {
		printf("FAILED: %s(%s)\n", cle.what(), ErrorString (cle.err()));
		return false;
	}

	printf ("done.\n    Builing program ... "); fflush (stdout);
	
	try {
		m_error = m_prg.build(m_devs);
		printf("done.\n");
	} catch (cl::Error cle) {
		printf("FAILED. Check logfile!\n");
		std::cout << "Build Status: "   << m_prg.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(m_devs[0])  << std::endl;
		std::cout << "Build Options:\t" << m_prg.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(m_devs[0]) << std::endl;
		std::cout << "Build Log:\t "    << m_prg.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_devs[0])     << std::endl;
		return false;
	}
	
	return true;
	
}



void
GPUSimulator::Simulate     () {

	if (!BuildProgram ("/usr/local/lib/GPUSimulator.cl")) 
		return;

	if (!PrepareKernel())
		return;

	SetDeviceData ();
	RunKernel ();
	GetDeviceData ();

}



bool
GPUSimulator::PrepareKernel () {
	
    printf ("    Forming kernel ... ");  fflush(stdout);
	
    try {
        m_kernel = cl::Kernel(m_prg, "Simulate", &m_error);
		printf ("done.\n"); fflush(stdout);
    } catch (cl::Error cle) {
        printf("FAILED: %s(%s)\n", cle.what(), ErrorString (cle.err()));
		return false;
    }

	return true;	

}


void 
GPUSimulator::SetDeviceData () {

	unsigned n[4] = {m_nt, m_nc, m_na, m_ne};

	
    printf("  Creating device arrays ... ");  fflush(stdout);
	ocl_tb1  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY, 2 * sizeof(float) *   m_sb->tb1->Size(), NULL, &m_error);
	ocl_sb1  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY, 2 * sizeof(float) *   m_sb->sb1->Size(), NULL, &m_error);
	ocl_agr  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->agr->Size(), NULL, &m_error);
	ocl_tr   = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *    m_sb->tr->Size(), NULL, &m_error);
	ocl_sr   = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *    m_sb->sr->Size(), NULL, &m_error);
	ocl_tb0  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->tb0->Size(), NULL, &m_error);
	ocl_sb0  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->sb0->Size(), NULL, &m_error);
	ocl_tmxy = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY, 2 * sizeof(float) *  m_sb->tmxy->Size(), NULL, &m_error);
	ocl_tmz  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->tmz->Size(), NULL, &m_error);
	ocl_smxy = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY, 2 * sizeof(float) *  m_sb->smxy->Size(), NULL, &m_error);
	ocl_smz  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->smz->Size(), NULL, &m_error);
	ocl_jac  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float) *   m_sb->jac->Size(), NULL, &m_error);
	ocl_gdt  = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY,     sizeof(float)                      , NULL, &m_error);
	ocl_n    = cl::Buffer (m_ctxt,  CL_MEM_READ_ONLY, 4 * sizeof(  int)                      , NULL, &m_error);
	ocl_rf   = cl::Buffer (m_ctxt, CL_MEM_WRITE_ONLY, 2 * sizeof(float) *    m_sb->rf->Size(), NULL, &m_error);
    ocl_mxy  = cl::Buffer (m_ctxt, CL_MEM_WRITE_ONLY, 2 * sizeof(float) *   m_sb->mxy->Size(), NULL, &m_error);
    ocl_mz   = cl::Buffer (m_ctxt, CL_MEM_WRITE_ONLY,     sizeof(float) *    m_sb->mz->Size(), NULL, &m_error);

	printf("done.\n  Pushing data to device ... ");  fflush(stdout);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_tb1, CL_TRUE, 0, 2 * sizeof(float) *  m_sb->tb1->Size(),  &m_sb->tb1->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_sb1, CL_TRUE, 0, 2 * sizeof(float) *  m_sb->sb1->Size(),  &m_sb->sb1->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_agr, CL_TRUE, 0,     sizeof(float) *  m_sb->agr->Size(),  &m_sb->agr->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (  ocl_tr, CL_TRUE, 0,     sizeof(float) *   m_sb->tr->Size(),   &m_sb->tr->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (  ocl_sr, CL_TRUE, 0,     sizeof(float) *   m_sb->sr->Size(),   &m_sb->sr->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_tb0, CL_TRUE, 0,     sizeof(float) *  m_sb->tb0->Size(),  &m_sb->tb0->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_sb0, CL_TRUE, 0,     sizeof(float) *  m_sb->sb0->Size(),  &m_sb->sb0->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (ocl_tmxy, CL_TRUE, 0, 2 * sizeof(float) * m_sb->tmxy->Size(), &m_sb->tmxy->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_tmz, CL_TRUE, 0,     sizeof(float) *  m_sb->tmz->Size(),  &m_sb->tmz->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (ocl_smxy, CL_TRUE, 0, 2 * sizeof(float) * m_sb->smxy->Size(), &m_sb->smxy->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_smz, CL_TRUE, 0,     sizeof(float) *  m_sb->smz->Size(),  &m_sb->smz->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_jac, CL_TRUE, 0,     sizeof(float) *  m_sb->jac->Size(),  &m_sb->jac->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_gdt, CL_TRUE, 0,     sizeof(float)                     ,             &m_gdt, NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (   ocl_n, CL_TRUE, 0, 4 * sizeof(int)                       ,              &n[0], NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (  ocl_rf, CL_TRUE, 0, 2 * sizeof(float) *   m_sb->rf->Size(),   &m_sb->rf->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer ( ocl_mxy, CL_TRUE, 0, 2 * sizeof(float) *  m_sb->mxy->Size(),  &m_sb->mxy->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueWriteBuffer (  ocl_mz, CL_TRUE, 0,     sizeof(float) *   m_sb->mz->Size(),   &m_sb->mz->At(0), NULL, &m_event);

	printf("done.\n  Assigning arguments to kernel ... ");  fflush(stdout);
	m_error = m_kernel.setArg( 0, ocl_tb1 );
	m_error = m_kernel.setArg( 1, ocl_sb1 );
	m_error = m_kernel.setArg( 2, ocl_agr );
	m_error = m_kernel.setArg( 3, ocl_tr  );
	m_error = m_kernel.setArg( 4, ocl_sr  );
	m_error = m_kernel.setArg( 5, ocl_tb0 );
	m_error = m_kernel.setArg( 6, ocl_sb0 );
	m_error = m_kernel.setArg( 7, ocl_tmxy);
	m_error = m_kernel.setArg( 8, ocl_tmz );
	m_error = m_kernel.setArg( 9, ocl_smxy);
	m_error = m_kernel.setArg(10, ocl_smz );
	m_error = m_kernel.setArg(11, ocl_jac );
	m_error = m_kernel.setArg(12, ocl_gdt );
	m_error = m_kernel.setArg(13, ocl_n   );
	m_error = m_kernel.setArg(14, ocl_rf  );
	m_error = m_kernel.setArg(15, ocl_mxy );
	m_error = m_kernel.setArg(16, ocl_mz  );

	printf ("done.\n");
    m_cmdq.finish();

}


void 
GPUSimulator::GetDeviceData () {

	printf("  Getting results back from kernel ... ");  fflush(stdout);
    m_error = m_cmdq.enqueueReadBuffer(ocl_rf, CL_TRUE, 0, 2 * sizeof(float) *  m_sb->rf->Size(),  &m_sb->rf->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueReadBuffer(ocl_rf, CL_TRUE, 0, 2 * sizeof(float) * m_sb->mxy->Size(), &m_sb->mxy->At(0), NULL, &m_event);
    m_error = m_cmdq.enqueueReadBuffer(ocl_rf, CL_TRUE, 0,     sizeof(float) *  m_sb->mz->Size(),  &m_sb->mz->At(0), NULL, &m_event);
	printf("done\n");  

}


void 
GPUSimulator::RunKernel () {

	printf("  Running kernel ... "); fflush (stdout);
	
    m_error = m_cmdq.enqueueNDRangeKernel (m_kernel, cl::NullRange, cl::NDRange(m_na), cl::NullRange, NULL, &m_event); 
    printf("FAILED: [clEnqueueNDRangeKernel: %s]", ErrorString(m_error));  fflush (stdout);
	
    m_cmdq.finish();
	
	printf ("... done.\n");
	
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

