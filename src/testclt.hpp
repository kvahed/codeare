#include "options.h"
#include "RemoteConnector.hpp"
#include "LocalConnector.hpp"
#include "modules/FFT.hpp"
#include "MatrixOperations.hpp"

#ifndef __WIN32__
    #include "config.h"
#endif

#ifndef SVN_REVISION
	#define SVN_REVISION "unkown"
#endif

using namespace std;
using namespace RRClient;

#include <time.h>
#include <stdio.h>

char*  name;
char*  base;
char*  data;
char*  config;
char*  verbose;
char*  test;
bool   remote;
bool   pulses;

#include "tests/tests.hpp"

bool init (int argc, char** argv);


