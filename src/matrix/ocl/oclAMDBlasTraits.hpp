# ifndef __OCL_AMD_BLAS_TRAITS_HPP__



  /************
   ** makros **
   ************/  
  # define __OCL_AMD_BLAS_TRAITS_HPP__



  /**************
   ** includes **
   **************/

  // AMD BLAS
  # include <clAmdBlas.h>
  
  # include "timer.h"




  /***********************
   ** Trait definitions **
   ***********************/
  
  
  /**********************************
   * Traits for AMD Blas operations *
   *  (base struct)                 *
   **********************************/  
  template <class T>
  struct amdBlasTraits
  { /* -- */ };
  
  
  
  /****************************************
   * AMD Blas operations single precision *
   *  (derived)                           *
   ****************************************/
  template <>
  struct amdBlasTraits <float>
  {
  
  
    typedef    float elem_type;
    typedef cl_float   cl_type;

    
    /**
     * @brief           getters for factors
     */
    static cl_type One  ()
    { return 1.; }
    static cl_type Zero ()
    { return 0.; }

  
    /**
     * @brief           General matrix-matrix product.
     */
    static
    cl_int
    GEMM                (clAmdBlasOrder order,
                         clAmdBlasTranspose transA, clAmdBlasTranspose transB,
                         size_t M, size_t N, size_t K,
                         cl_type alpha, const cl_mem A, size_t lda,
                                        const cl_mem B, size_t ldb,
                         cl_type  beta,       cl_mem C, size_t ldc,
                         cl_uint numCommandQueues, cl_command_queue * commandQueues,
                         cl_uint numEventsInWaitList, const cl_event * eventWaitList, cl_event * events)
    {
      
      clAmdBlasSgemm (order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, numCommandQueues, commandQueues,
                      numEventsInWaitList, eventWaitList, events);
      
    }
    
  
  };
  
  
  /****************************************
   * AMD Blas operations double precision *
   *  (derived)                           *
   ****************************************/
  template <>
  struct amdBlasTraits <double>
  {
  
  
    typedef    double elem_type;
    typedef cl_double   cl_type;
    
    
    /**
     * @brief           getters for factors
     */
    static cl_type One  ()
    { return 1.; }
    static cl_type Zero ()
    { return 0.; }
  
    
    /**
     * @brief           General matrix-matrix product.
     */
    static
    cl_int
    GEMM                (clAmdBlasOrder order,
                         clAmdBlasTranspose transA, clAmdBlasTranspose transB,
                         size_t M, size_t N, size_t K,
                         cl_type alpha, const cl_mem A, size_t lda,
                                         const cl_mem B, size_t ldb,
                         cl_type  beta,       cl_mem C, size_t ldc,
                         cl_uint numCommandQueues, cl_command_queue * commandQueues,
                         cl_uint numEventsInWaitList, const cl_event * eventWaitList, cl_event * events)
    {
      
      
      // TODO: use EXtended version
      clAmdBlasDgemm (order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, numCommandQueues, commandQueues,
                      numEventsInWaitList, eventWaitList, events);
      
    }
    
  
  };
  
  
  
  /************************************************
   * AMD Blas operations complex single precision *
   *  (derived)                                   *
   ************************************************/
  template <>
  struct amdBlasTraits <cxfl>
  {
  
  
    typedef      cxfl elem_type;
    typedef cl_float2   cl_type;
    
    
    /**
     * @brief           getters for factors
     */
    static cl_type One  ()
    { cl_type one = {1., 0.}; return one; }
    static cl_type Zero ()
    { cl_type zero = {0., 0.}; return zero; }
  
    
    /**
     * @brief           General matrix-matrix product.
     */
    static
    cl_int
    GEMM                (clAmdBlasOrder order,
                         clAmdBlasTranspose transA, clAmdBlasTranspose transB,
                         size_t M, size_t N, size_t K,
                         cl_type alpha, const cl_mem A, size_t lda,
                                         const cl_mem B, size_t ldb,
                         cl_type  beta,       cl_mem C, size_t ldc,
                         cl_uint numCommandQueues, cl_command_queue * commandQueues,
                         cl_uint numEventsInWaitList, const cl_event * eventWaitList, cl_event * events)
    {
      
      
      // TODO: use EXtended version
      clAmdBlasCgemm (order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, numCommandQueues, commandQueues,
                      numEventsInWaitList, eventWaitList, events);
      
    }
    
  
  };
  


  /************************************************
   * AMD Blas operations complex double precision *
   *  (derived)                                   *
   ************************************************/
  template <>
  struct amdBlasTraits <cxdb>
  {
  
  
    typedef       cxdb elem_type;
    typedef cl_double2   cl_type;
    
    
    /**
     * @brief           getters for factors
     */
    static cl_type One  ()
    { cl_type one = {1., 0.}; return one; }
    static cl_type Zero ()
    { cl_type zero = {0., 0.}; return zero; }
  
    
    /**
     * @brief           General matrix-matrix product.
     */
    static
    cl_int
    GEMM                (clAmdBlasOrder order,
                         clAmdBlasTranspose transA, clAmdBlasTranspose transB,
                         size_t M, size_t N, size_t K,
                         cl_type alpha, const cl_mem A, size_t lda,
                                         const cl_mem B, size_t ldb,
                         cl_type  beta,       cl_mem C, size_t ldc,
                         cl_uint numCommandQueues, cl_command_queue * commandQueues,
                         cl_uint numEventsInWaitList, const cl_event * eventWaitList, cl_event * events)
    {
      
      
      // TODO: use EXtended version
      clAmdBlasZgemm (order, transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc, numCommandQueues, commandQueues,
                      numEventsInWaitList, eventWaitList, events);
      
    }
    
  
  };




# endif // __OCL_AMD_BLAS_TRAITS_HPP__
