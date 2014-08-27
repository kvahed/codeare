/*
 * Assert.hpp
 *
 *  Created on: Aug 20, 2014
 *      Author: kvahed
 */

#include <assert.h>


#ifndef __ASSERT_HPP__
#define __ASSERT_HPP__


#define stringize(s) #s
#define XSTR(s) stringize(s)
#if !defined NDEBUG
#pragma warning (disable : 4273)
	void abort (void);
#pragma warning (default : 4273)
# if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
# define op_assert(a, b, c) \
do { \
    if (0 == (a)) { \
    	fprintf(stderr, "    ERROR: assertion failed: %s, %s(), %d at \'%s\'\n", \
		    __FILE__, __func__, __LINE__, XSTR(a)); \
		    std::cerr << "           Matrix dimensions do not match: (" << \
		    		size(b) << ") != (" << size(c) << ")\n"; \
        abort(); \
    } \
} while (0)
# else
# define op_assert(a, b, c) \
do { \
	if (0 == (a)) { \
		fprintf(stderr, "    ERROR: assertion failed: %s, %d at \'%s\'\n", \
			__FILE__, __LINE__, XSTR(a)); \
		std::cerr << "           Matrix dimensions do not match: (" << \
				size(b) << ") != (" << size(c) << ")\n"; \
		abort(); \
	} \
} while (0)
# endif
#else
# define op_assert(a, b, c) (void)0
#endif



#endif /* __ASSERT_HPP__ */
