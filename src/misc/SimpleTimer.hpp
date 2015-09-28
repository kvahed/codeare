#ifndef __SIMPLE_TIMER_HPP__
#define __SIMPLE_TIMER_HPP__

#ifdef _MSC_VER 
#include <Windows.h>
#else         
#include <sys/time.h>
#endif

#include <iostream>
#include <boost/format.hpp>

#include <boost/timer/timer.hpp>
#include <cmath>


class SimpleTimer {

public:

    inline SimpleTimer () : _verbose(false) { 
        _timer.start();
    }

    inline SimpleTimer (const std::string& identifier,
                        std::ostream& os = std::cout) : _os(&os), _verbose(true) { 
        *_os << "Processing " << identifier.c_str() << " ... \n";
        _timer.start();
    }

    inline ~SimpleTimer () {
        if (!_timer.is_stopped()) 
            _timer.stop();
        if (_verbose)
            *_os << "... done. WTime: " << boost::format("%.4f") % (_timer.format()) << "\n\n";
    }

    inline void Stop () {
        _timer.stop();
    }

    inline void Resume () {
        _timer.resume();
    }

    inline boost::timer::cpu_times Elapsed () const {
        return _timer.elapsed();
    }

    inline std::string Format () const {
        return _timer.format();
    }

private:

    boost::timer::cpu_timer _timer;
    std::ostream  *_os;
    bool _verbose;
    
};


/*
class SimpleTimer {

#ifdef _MSC_VER
    typedef LARGE_INTEGER TimeType;
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

*/

#endif // __SIMPLE_TIMER_HPP__
