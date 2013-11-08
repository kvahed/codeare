
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
  
//    if (get_work_dim () == 1)
//    {
      *index = get_local_id (0);
//    }
//    else
//    {
//      return false;
//    }
   
    *global_size = get_global_size (0);
    
    *global_inc = *global_size;
    
    *local_position = *index;
  
    return true;
  
  }

# endif // __OCL_KERNEL_HEADER__




/**
 * @brief                 Elementwise increment vector.
 *
 * @param  arg1           Start address of first vector.
 * @param  inc            Scalar value (increment).
 * @param  num_elems      Number of elements.
 */
__kernel
void
assign                    ( __global A_type_n *      arg1,
                            __global A_type   *    scalar,
                            __global      int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems / vec_len; i += global_inc)
  {
    arg1 [i] = *scalar;
  }
 
}



/**
 * @brief                 Bitwise AND (mask).
 *
 * @param  arg1           Vector to be masked.
 * @param  mask           Complex mask vector.
 * @param  res            Cross_section or zero.
 * @param  num_elems      Number of vectors' elements.
 */
__kernel
void
bitw_and                  ( __global A_type *      arg1,
                            __global   bool *      mask,
                            __global A_type *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (mask [i] ? arg1 [i] : res [i]);
  }

}



/**
 * @brief                 Scalar comparison.
 *
 * @param  arg1           Matrix.
 * @param  scalar         Scalar.
 * @param  res            Result matrix.
 * @param  num_elems      Number of elements.
 */
__kernel
void
scalar_equal              ( __global A_type *      arg1,
                            __global A_type *    scalar,
                            __global   bool *       res,
                            __global    int * num_elems )
{
  
  __private A_type s = *scalar;

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (arg1 [i].x == s.x) && (arg1 [i].y == s.y);
  }
  
}



/**
 * @brief                 Scalar comparison.
 *
 * @param  arg1           Matrix.
 * @param  scalar         Scalar.
 * @param  res            Result matrix.
 * @param  num_elems      Number of elements.
 */
__kernel
void
scalar_inequal              ( __global A_type *      arg1,
                              __global A_type *    scalar,
                              __global   bool *       res,
                              __global    int * num_elems )
{
  
  __private A_type s = *scalar;

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (arg1 [i].x != s.x) || (arg1 [i].y != s.y);
  }
  
}



/**
 * @brief                 Elementwise comparison.
 *
 * @param  arg1           1st Matrix.
 * @param  arg2           2nd Matrix.
 * @param  res            Result matrix.
 * @param  num_elems      Number of elements.
 */
__kernel
void
vector_equal              ( __global A_type *      arg1,
                            __global A_type *      arg2,
                            __global   bool *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (arg1 [i].x == arg2 [i].x) && (arg1 [i].y == arg2 [i].y);
  }
  
}



/**
 * @brief                 Elementwise comparison.
 *
 * @param  arg1           1st Matrix.
 * @param  arg2           2nd Matrix.
 * @param  res            Result matrix.
 * @param  num_elems      Number of elements.
 */
__kernel
void
vector_inequal            ( __global A_type *      arg1,
                            __global A_type *      arg2,
                            __global   bool *       res,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    res [i] = (arg1 [i].x != arg2 [i].x) || (arg1 [i].y != arg2 [i].y);
  }
  
}



/**
 * @brief                 Deep copy of a vector.
 *
 * @param  arg1           Destination.
 * @param  arg2           Source.
 * @param  num_elems      Number of elements.
 */
__kernel
void
copy_buffer               ( __global A_type *      arg1,
                            __global A_type *      arg2,
                            __global    int * num_elems )
{

  int index, global_size, global_inc, local_position;
  
  /* initialize parameters */
  init_params (&index, &global_size, &global_inc, &local_position);

  /* calculation */
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg1 [i] = arg2 [i];
  }
  
}
