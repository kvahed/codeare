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



__kernel void add (__global float * arg1, __global float * arg2, __global float * result, __global int * size)
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
  
  for (int i = local_position; i < *size; i += global_inc)
  {
    result [i] = arg1 [i] + arg2 [i];
  }
 
}


typedef float elem_type;

#pragma OPENCL_EXTENSION cl_amd_printf : enable
__kernel void convolution1 (__global elem_type * signal, __global int * singal_length_x, __global int * signal_length_y,
                            __global elem_type * filter, __global int * filter_length_x, __global int * filter_length_y)
{

  int2 index;

  // identify thread position
  if (get_work_dim() == 2)
  {
    index.x = get_global_id (0), index.y = get_global_id (1);
  }
  else
  {
    return;
  }

  // load signal block to shared memory
  
  // execute convolution
  
  // write back calculated pixel value

}
