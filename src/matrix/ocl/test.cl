__kernel void testKernel (__global float * numbers, __global float * result) {

  float tmp_number = min (*numbers, *(numbers + 1));

  *numbers = max (*numbers, *(numbers + 1));

  *(numbers + 1) = tmp_number;

}


__kernel void testKernel2 (__global float * numbers, __global float * result) {

  float tmp_number = max (*numbers, *(numbers + 1));

  *numbers = min (*numbers, *(numbers + 1));

  *(numbers + 1) = tmp_number;

}


/**
 * @brief                 Elementwsie add elements of the two vectors.
 *
 * @param  arg1           Start address of first vector.
 * @param  arg2           Start address of second vector.
 * @param  result         Start address of result vector.
 * @param  num_elems      Number of elements of all vectors.
 */
__kernel
void
add                       ( __global float *      arg1,
                            __global float *      arg2,
                            __global float *    result,
                            __global   int  * num_elems )
{

  int index;

  if (get_work_dim () == 1)
  {
    index = get_local_id (0);
  }
  else
  {
    return;
  }
 
  int global_size = get_global_size (0);
  
  int global_inc = global_size;
  
  int local_position = index;
  
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    result [i] = arg1 [i] + arg2 [i];
  }
  
/*  result [0] = 4.44;
  result [(*num_elems)/2] = 8.88;
  result [(*num_elems)] = 16.1616;
*/ 
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
inc                       ( __global float *      arg1,
                            __global float *       inc,
                            __global   int * num_elems )
{

  int index;
  
  __local float scalar;
  scalar = *inc;

  if (get_work_dim () == 1)
  {
    index = get_local_id (0);
  }
  else
  {
    return;
  }
 
  int global_size = get_global_size (0);
  
  int global_inc = global_size;
  
  int local_position = index;
  
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg1 [i] = arg1 [i] + scalar;
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
assign                    ( __global float *      arg1,
                            __global float *    scalar,
                            __global   int * num_elems )
{

  int index;
  
  if (get_work_dim () == 1)
  {
    index = get_local_id (0);
  }
  else
  {
    return;
  }
 
  int global_size = get_global_size (0);
  
  int global_inc = global_size;
  
  int local_position = index;
  
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    arg1 [i] = *scalar;
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
copy_buffer               ( __global float *      arg1,
                            __global float *      arg2,
                            __global   int * num_elems )
{

  int index;

  if (get_work_dim () == 1)
  {
    index = get_local_id (0);
  }
  else
  {
    return;
  }
 
  int global_size = get_global_size (0);
  
  int global_inc = global_size;
  
  int local_position = index;
  
  for (int i = local_position; i < *num_elems; i += global_inc)
  {
    
    arg1 [i] = arg2 [i];

  }
  
}



//typedef float elem_type;

//#pragma OPENCL_EXTENSION cl_amd_printf : enable
//__kernel void convolution1 (__global elem_type * signal, __global int * singal_length_x, __global int * signal_length_y,
//                            __global elem_type * filter, __global int * filter_length_x, __global int * filter_length_y)
//{

//  int2 index;

  // identify thread position
//  if (get_work_dim() == 2)
//  {
//    index.x = get_global_id (0), index.y = get_global_id (1);
//  }
  //else
  //{
  //  return;
  //}

  // load signal block to shared memory
  
  // execute convolution
  
  // write back calculated pixel value

//}
