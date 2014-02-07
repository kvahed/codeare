#ifndef __SIMPLE_TIMER_HPP__
#define __SIMPLE_TIMER_HPP__

#ifdef _MSC_VER 
#include <windows.h>
#else         
#include <sys/time.h>
#endif

#include <iostream>
#include <boost/format.hpp>

class SimpleTimer {

#ifdef _MSC_VER
    typedef LONG_INTEGER TimeType ;
#else
    typedef timeval      TimeType;
#endif
    
public: 
    inline SimpleTimer (const std::string& identifier = "", std::ostream& os = std::cout) :
        _os(&os), _stopped(false), _freq (1.0e-6), _net(0.) {
        *_os << "Processing " << identifier.c_str() << " ... \n";
        Resume();
    }
    
    inline ~SimpleTimer() {
        if (!_stopped)
            Stop();
        *_os << "... done. WTime: " << boost::format("%.4f") % (_net * _freq) << " ... \n\n";
    }

    inline void Stop () {
        _stopped = true; 
        TimeType tmp;
#ifdef _MSC_VER
        QueryPerformanceCounter(&tmp);
        _net += tmp.QuadPart - _start.QuadPart;
#else
        gettimeofday(&tmp, NULL);
        _net += (1.0e6 * tmp.tv_sec  + tmp.tv_usec) - (1.0e6 * _start.tv_sec  + _start.tv_usec);
#endif
    }

    inline void Resume () {
        _stopped = false;
#ifdef _MSC_VER
        QueryPerformanceCounter(&_start);
#else
        gettimeofday(&_start, NULL);
#endif
    }
    
private: 

    bool      _stopped;
    std::ostream  *_os;
    TimeType  _start;
    double    _net;
    double    _freq;
    
};

#endif // __SIMPLE_TIMER_HPP__
