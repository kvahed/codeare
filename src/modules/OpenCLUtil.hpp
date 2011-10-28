//from http://www.songho.ca/opengl/gl_vbo.html

#ifndef __OPENCL_UTIL_H__
#define __OPENCL_UTIL_H__

char*
file_contents  (const char *filename, int *length);

const char* 
oclErrorString (cl_int error);

#endif // __OPENCL_UTIL_H__

