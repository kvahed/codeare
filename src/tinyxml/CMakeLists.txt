if (!"${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()
add_library (tinyxml tinystr.cpp tinystr.h tinyxml.cpp tinyxml.h tinyxmlerror.cpp tinyxmlparser.cpp xpath_processor.h xpath_processor.cpp xpath_expression.cpp xpath_expression.h xpath_static.cpp xpath_static.h action_store.cpp action_store.h xpath_stack.cpp xpath_stack.h lex_util.cpp lex_util.h xpath_stream.cpp xpath_stream.h node_set.cpp node_set.h xml_util.cpp xml_util.h xpath_syntax.cpp xpath_syntax.h tokenlist.cpp tokenlist.h)
