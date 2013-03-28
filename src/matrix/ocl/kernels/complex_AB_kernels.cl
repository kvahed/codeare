
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
    arg1 [i] = arg1 [i] + scalar;
  }
 
}



/**
 * @brief                 Elementwise multiplication with complex scalar.
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

  // storage for real part
  __private A_type re;

  // calculation
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    re = arg1 [i] * *scalar;
    res [i].x = re.x - re.y;
    res [i].y = dot (arg1 [i].yx, *scalar);
  }
  
}



/**
 * @brief                 Elementwise multiplication with complex vector.
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

  // storage for real part
  __private A_type re;

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    re = arg1 [i] * arg2 [i];
    res [i].x = re.x - re.y;
    res [i].y = dot (arg1 [i].yx, arg2 [i]);
  }
  
}



/**
 * @brief                 Elementwise division by complex scalar.
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

  __private A_type re;
  __private A_type div;

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    re = arg1 [i].yx * *scalar;
    div.x = pow (length (*scalar), 2);
    if (div.x == 0)
    {
      res [i].x = 0;
      res [i].y = 0;
    }
    else
    {
      res [i].x = (dot (arg1 [i], *scalar)) / div.x;
      res [i].y = (re.x - re.y) / div.x;
    }  
  }
  
}



/**
 * @brief                 Elementwise division by complex vector.
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

  __private A_type re;
  __private A_type div;

  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    re = arg1 [i].yx * arg2 [i];
    div.x = pow (length (arg2 [i]), 2);
    if (div.x == 0)
    {
      res [i].x = 0;
      res [i].y = 0;
    }
    else
    {
      res [i].x = (dot (arg1 [i], arg2 [i])) / div.x;
      res [i].y = (re.x - re.y) / div.x;
    }  
  }
  
}



/**
 * @brief                 Complex Matrix-Matrix product (not inplace).
 *
 * @param  arg1           First matrix.
 * @param  arg2           Second matrix.
 * @param  res            Result matrix.
 * @param  m              Number of rows of first matrix.
 * @param  n              Number of cols of first matrix
 *                                  rows of second matrix.
 * @param  k              Number of cols of second matrix.
 * @param  transA         Transpose first matrix if inequals 0.
 * @param  transB         Transpose second matrix if inequals 0.
 */
__kernel
void
mat_prod                  ( __global A_type *   arg1,
                            __global B_type *   arg2,
                            __global A_type *    res,
                            __global    int *      m,
                            __global    int *      n,
                            __global    int *      k,
                            __global    int * transA,
                            __global    int * transB )
{

  

}



/**
 * @brief                 "Elementwise" cast of vector.
 *
 * @param  arg1           Vector to be casted.
 * @param  arg2           Vector containing cast results.
 * @param  num_elems      Number of elements.
 */
/*__kernel
void
vector_cast               ( __global A_type *      arg1,
                            __global B_type *      arg2,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  // initialize parameters //
  init_params (&index, &global_size, &global_inc, &local_position);
  
  // calculation //
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg2 [i] = (B_type) arg1 [i];
  }

}*/
