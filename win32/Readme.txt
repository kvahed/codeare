* Make directory C:\codeare

* From http://sourceforge.net/projects/omniorb/files/omniORB download
  Version 4.1.1 for Windows vc6. (IMPORTANT! Need perfekt match.)
* Unzip in C:\codeare 
* Copy all files in omniORB\lib\x86_win32 Z:\n4\x86\extsw\lib

* From http://www.microsoft.com/en-us/download/details.aspx?id=29
  download Microsoft Visual C++ 2008 Redistributable Package (x86)

* From http://slproweb.com/products/Win32OpenSSL.html download OpenSSL
  for developers http://slproweb.com/download/Win32OpenSSL-1_0_1e.exe
* Install

* Download current boost and unpack somewhere
* VC6 bug prevents compilation so one file needs one tiny edit:
  In boost/assert.hpp replace "std::abort()" by "exit(-1)"

