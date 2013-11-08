
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
