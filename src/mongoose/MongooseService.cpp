/*
 * MongooseService.cpp
 *
 *  Created on: Aug 7, 2013
 *      Author: kvahed
 */

#include "MongooseService.hpp"
#include "Workspace.hpp"

#include <boost/thread.hpp>
#include <boost/bind.hpp>

namespace codeare {
namespace service {


MongooseService* MongooseService::_instance = 0;


MongooseService& MongooseService::Instance () {

	if (_instance == 0)
		_instance = new MongooseService();

	return *_instance;

}

static const std::string head = "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\nContent-Length: %d\r\n\r\n%s";

static int handler (struct mg_connection *conn) {

	mg_get_request_info(conn);
	std::stringstream wd;
	wd << wspace;
	std::string ws = wd.str();

	mg_printf (conn, head.c_str(), ws.length(), ws.c_str());

	return ws.length();

}

static void Serve () {

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
    //boost::thread bt (Serve);
}

MongooseService::~MongooseService() {}


} /* namespace service */
} /* namespace codeare */
