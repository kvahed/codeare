
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
 * @brief                 Elementwise add elements of the two vectors (complex, real).
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
    result [i].x = arg1 [i].x + arg2 [i];
    result [i].y = arg1 [i].y;
  }

}



/**
 * @brief                 Elementwise subtract elements of the two vectors (complex, real).
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
    result [i].x = arg1 [i].x - arg2 [i];
    result [i].y = arg1 [i].y;
  }

}



/**
 * @brief                 Elementwise real increment complex vector.
 *
 * @param  arg1           Start address of first vector.
 * @param  inc            Scalar value (increment).
 * @param  num_elems      Number of elements.
 */
__kernel
void
inc                       ( __global A_type   *      arg1,
                            __global B_type   *       inc,
                            __global      int * num_elems )
{
  
  __local B_type scalar;
  scalar = *inc;

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg1 [i].x = arg1 [i].x + scalar;
  }
 
}



/**
 * @brief                 Elementwise multiplication with real scalar.
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
  
  // initialize parameters
  init_params (&index, &global_size, &global_inc, &local_position);

  // calculation
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i].x = arg1 [i].x * *scalar;
    res [i].y = arg1 [i].y * *scalar;
  }
  
}



/**
 * @brief                 Elementwise multiplication with real vector.
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
  
  // initialize parameters //
  init_params (&index, &global_size, &global_inc, &local_position);

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i].x = arg1 [i].x * arg2 [i];
    res [i].y = arg1 [i].y * arg2 [i];
  }
  
}



/**
 * @brief                 Elementwise division by real scalar.
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
  
  // initialize parameters //
  init_params (&index, &global_size, &global_inc, &local_position);

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    if (*scalar == 0)
    {
      res [i].x = 0;
      res [i].y = 0;
    }
    else
    {
      res [i].x = arg1 [i].x / *scalar;
      res [i].y = arg1 [i].y / *scalar;
    }  
  }
  
}



/**
 * @brief                 Elementwise division by real vector.
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
  
  // initialize parameters //
  init_params (&index, &global_size, &global_inc, &local_position);

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    if (arg2 [i] == 0)
    {
      res [i].x = 0;
      res [i].y = 0;
    }
    else
    {
      res [i].x = arg1 [i].x / arg2 [i];
      res [i].y = arg1 [i].y / arg2 [i];
    }  
  }
  
}



/**
 * @brief                 Elementwise assignment of complex vector by real scalar.
 *
 * @param  arg1           Start address of first vector.
 * @param  inc            Scalar value (increment).
 * @param  num_elems      Number of elements.
 */
__kernel
void
assign                    ( __global A_type   *      arg1,
                            __global B_type   *    scalar,
                            __global    int   * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems / vec_len; i += global_inc)
  {
    arg1 [i].x = *scalar;
    arg1 [i].y = 0;
  }
 
}

