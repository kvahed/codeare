#ifndef __SIMPLE_TIMER_HPP__
#define __SIMPLE_TIMER_HPP__

#include "Toolbox.hpp"
#include "cycle.h"            // FFTW cycle implementation

#include <iostream>
#include <boost/format.hpp>


class SimpleTimer {

public: 
    inline SimpleTimer (const std::string& identifier = "", std::ostream& os = std::cout) :
        _start (getticks()), _identifier(""), _os(&os), _stopped(false) , _net(0) {
        *_os << "Processing " << identifier.c_str() << " ... \n";
    }

    inline ~SimpleTimer() {
        if (!_stopped)
            Stop();
        float e = 1.0 / Toolbox::Instance()->ClockRate() * _net;
        *_os << " ... done. WTime: " << boost::format("%.4f") % e << " ... \n\n";
    }

    inline void Stop () {
        _net += elapsed(getticks(), _start);
        _stopped = true;
    }

    inline void Resume () {
        _stopped = false;
        _start = getticks();
    }
    
    const std::string _identifier;

private: 
    ticks _start;
    ticks _net;
    bool _stopped;
    ostream *_os;
    
};

#endif // __SIMPLE_TIMER_HPP__
