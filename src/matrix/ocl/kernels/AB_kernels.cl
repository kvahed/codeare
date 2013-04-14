
/**
 * added by oclConnection:
 *    A_type, A_type_n
 *    B_type, B_type_n
 *    vec_len
 */

# ifndef __OCL_KERNEL_HEADER__

  # define __OCL_KERNEL_HEADER__
  
  /**
   * @brief                 Helper function for initializing kernel parameters.
   *
   * @return                State of success.
   */
  bool
  init_params               ( int * index,
                              int * global_size,
                              int * global_inc,
                              int * local_position )
  {
  
    if (get_work_dim () == 1)
    {
      *index = get_local_id (0);
    }
    else
    {
      return false;
    }
   
    *global_size = get_global_size (0);
    
    *global_inc = *global_size;
    
    *local_position = *index;
  
    return true;
  
  }

# endif // __OCL_KERNEL_HEADER__



/**
 * @brief                 Elementwise add elements of the two vectors.
 *
 * @param  arg1           Start address of first vector.
 * @param  arg2           Start address of second vector.
 * @param  result         Start address of result vector.
 * @param  num_elems      Number of elements of all vectors.
 */
__kernel
void
vector_add                ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global A_type *    result,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    result [i] = arg1 [i] + arg2 [i];
  }

}



/**
 * @brief                 Elementwise subtract elements of the two vectors.
 *
 * @param  arg1           Start address of first vector.
 * @param  arg2           Start address of second vector.
 * @param  result         Start address of result vector.
 * @param  num_elems      Number of elements of all vectors.
 */
__kernel
void
vector_sub                ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global A_type *    result,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    result [i] = arg1 [i] - arg2 [i];
  }

}



/**
 * @brief                 Elementwise increment vector.
 *
 * @param  arg1           Start address of first vector.
 * @param  inc            Scalar value (increment).
 * @param  num_elems      Number of elements.
 */
__kernel
void
inc                       ( __global A_type_n *      arg1,
                            __global B_type   *       inc,
                            __global      int * num_elems )
{
  
  __local B_type scalar;
  scalar = *inc;

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems / vec_len; i += global_inc)
  {
    arg1 [i] = arg1 [i] + scalar;
  }
 
}



/**
 * @brief                 Elementwise multiplication with scalar.
 *
 * @param  arg1           Vector.
 * @param  scalar         Factor scalar.
 * @param  res            Result vector.
 * @param  num_elems      Number of elements.
 */
__kernel
void
scalar_mult               ( __global A_type *      arg1,
                            __global B_type *    scalar,
                            __global A_type *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = arg1 [i] * *scalar;
  }
  
}



/**
 * @brief                 Elementwise multiplication with vector.
 *
 * @param  arg1           First vector.
 * @param  arg2           Second vector.
 * @param  res            Result vector.
 * @param  num_elems      Number of elements.
 */
__kernel
void
vector_mult               ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global A_type *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = ((A_type) arg2 [i]) * arg1 [i];
  }
  
}



/**
 * @brief                 Elementwise division by scalar.
 *
 * @param  arg1           Vector.
 * @param  scalar         Divisor scalar.
 * @param  res            Result vector.
 * @param  num_elems      Number of elements.
 */
__kernel
void
scalar_div                ( __global A_type *      arg1,
                            __global B_type *    scalar,
                            __global A_type *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (*scalar != 0) ? arg1 [i] / *scalar : 0;
  }
  
}



/**
 * @brief                 Elementwise division by vector.
 *
 * @param  arg1           First vector.
 * @param  arg2           Second vector.
 * @param  res            Result vector.
 * @param  num_elems      Number of elements.
 */
__kernel
void
vector_div                ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global A_type *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    A_type tmp = (A_type) arg2 [i];
    res [i] = (arg2 [i] != 0) ? arg1 [i] / tmp : 0;
  }
  
}



/**
 * @brief                 "Elementwise" cast of vector.
 *
 * @param  arg1           Vector to be casted.
 * @param  arg2           Vector containing cast results.
 * @param  num_elems      Number of elements.
 */
__kernel
void
vector_cast               ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);
  
  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg2 [i] = (B_type) arg1 [i];
  }

}
