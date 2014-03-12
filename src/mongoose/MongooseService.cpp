/*
 * MongooseService.cpp
 *
 *  Created on: Aug 7, 2013
 *      Author: kvahed
 */

#include "MongooseService.hpp"
#include "Workspace.hpp"

#ifdef HAVE_CXX11_THREAD
#  include <thread>
#  include <atomic>         // std::atomic
#  define thread_i std::thread
#else
#  include <boost/thread.hpp>
#  include <boost/bind.hpp>
#  define thread_i boost::thread
#endif
namespace codeare {
namespace service {


    MongooseService* MongooseService::_instance = 0;

    
    MongooseService& MongooseService::Instance () {
        
        if (_instance == 0)
            _instance = new MongooseService();
        
        return *_instance;
        
    }
    
    static const std::string head = "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\nContent-Length: %d\r\n\r\n%s";
    
    int handler (struct mg_connection *conn) {
        
        mg_get_request_info(conn);
        std::stringstream wd;
        wd << wspace;
        std::string ws = wd.str();
        
        mg_printf (conn, head.c_str(), ws.length(), ws.c_str());
        
        return ws.length();
        
    }
    
    void Serve () {
        struct mg_context *ctx;
        struct mg_callbacks callbacks;
        const char *options[] = {"listening_ports", wspace.p.Get<char*>("http_port"), NULL};
        
        memset(&callbacks, 0, sizeof(callbacks));
        callbacks.begin_request = handler;
        ctx = mg_start(&callbacks, NULL, options);
        getchar();
        mg_stop(ctx);
    }

    MongooseService::MongooseService() {
        thread_i thrd (Serve);
    }
    
    MongooseService::~MongooseService() {}
    
} /* namespace service */
} /* namespace codeare */
