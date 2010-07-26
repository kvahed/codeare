#ifndef __DLL_EXPORT_H__
#define __DLL_EXPORT_H__

#if defined (_MSC_VER)
    #define DLLEXPORT __declspec( dllexport )
#elif defined (__GNUC__)
    #define DLLEXPORT
#endif

#endif // __DLL_EXPORT_H__
