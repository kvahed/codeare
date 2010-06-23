
/*
 * Example daemon shell code for all of the requirements of a basic
 * linux daemon written in C.
 *
 * To use this code, search for 'TODO' and follow the directions.
 * 
 * To compile this file:
 *      gcc -o [daemonname] thisfile.c
 *
 * Substitute gcc with cc on some platforms.
 *
 * Peter Lombardo (peter AT lombardo DOT info)
 * 5/1/2006
 *
 */
 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <syslog.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
 
// TODO: Change '[daemonname]' to the name of _your_ daemon
#define DAEMON_NAME "jrrsd"
#define PID_FILE "/var/run/jrrsd.pid"
 
/**************************************************************************
    Function: Print Usage
 
    Description:
        Output the command-line options for this daemon.
 
    Params:
        @argc - Standard argument count
        @argv - Standard argument array
 
    Returns:
        returns void always
**************************************************************************/
void PrintUsage(int argc, char *argv[]) {
    if (argc >=1) {
        printf("Usage: %s -h -nn", argv[0]);
        printf("  Options:n");
        printf("      -ntDon't fork off as a daemon.n");
        printf("      -htShow this help screen.n");
        printf("n");
    }
}
 
/**************************************************************************
    Function: signal_handler
 
    Description:
        This function handles select signals that the daemon may
        receive.  This gives the daemon a chance to properly shut
        down in emergency situations.  This function is installed
        as a signal handler in the 'main()' function.
 
    Params:
        @sig - The signal received
 
    Returns:
        returns void always
**************************************************************************/
void signal_handler(int sig) {
 
    switch(sig) {
        case SIGHUP:
            syslog(LOG_WARNING, "Received SIGHUP signal.");
            break;
        case SIGTERM:
            syslog(LOG_WARNING, "Received SIGTERM signal.");
            break;
        default:
            syslog(LOG_WARNING, "Unhandled signal (%d) %s", strsignal(sig));
            break;
    }
}
 
/**************************************************************************
    Function: main
 
    Description:
        The c standard 'main' entry point function.
 
    Params:
        @argc - count of command line arguments given on command line
        @argv - array of arguments given on command line
 
    Returns:
        returns integer which is passed back to the parent process
**************************************************************************/
int main(int argc, char *argv[]) {
 
#if defined(DEBUG)
    int daemonize = 0;
#else
    int daemonize = 1;
#endif
 
    // Setup signal handling before we start
    signal(SIGHUP, signal_handler);
    signal(SIGTERM, signal_handler);
    signal(SIGINT, signal_handler);
    signal(SIGQUIT, signal_handler);
 
    int c;
    while( (c = getopt(argc, argv, "nh|help")) != -1) {
        switch(c){
            case 'h':
                PrintUsage(argc, argv);
                exit(0);
                break;
            case 'n':
                daemonize = 0;
                break;
            default:
                PrintUsage(argc, argv);
                exit(0);
                break;
        }
    }
 
    syslog(LOG_INFO, "%s daemon starting up", DAEMON_NAME);
 
    // Setup syslog logging - see SETLOGMASK(3)
#if defined(DEBUG)
    setlogmask(LOG_UPTO(LOG_DEBUG));
    openlog(DAEMON_NAME, LOG_CONS | LOG_NDELAY | LOG_PERROR | LOG_PID, LOG_USER);
#else
    setlogmask(LOG_UPTO(LOG_INFO));
    openlog(DAEMON_NAME, LOG_CONS, LOG_USER);
#endif
 
    /* Our process ID and Session ID */
    pid_t pid, sid;
 
    if (daemonize) {
        syslog(LOG_INFO, "starting the daemonizing process");
 
        /* Fork off the parent process */
        pid = fork();
        if (pid < 0) {
            exit(EXIT_FAILURE);
        }
        /* If we got a good PID, then
           we can exit the parent process. */
        if (pid > 0) {
            exit(EXIT_SUCCESS);
        }
 
        /* Change the file mode mask */
        umask(0);
 
        /* Create a new SID for the child process */
        sid = setsid();
        if (sid < 0) {
            /* Log the failure */
            exit(EXIT_FAILURE);
        }
 
        /* Change the current working directory */
        if ((chdir("/")) < 0) {
            /* Log the failure */
            exit(EXIT_FAILURE);
        }
 
        /* Close out the standard file descriptors */
        close(STDIN_FILENO);
        close(STDOUT_FILENO);
        close(STDERR_FILENO);
    }
 
    //****************************************************
    // TODO: Insert core of your daemon processing here
    //****************************************************
 
    syslog(LOG_INFO, "%s daemon exiting", DAEMON_NAME);
 
    //****************************************************
    // TODO: Free any allocated resources before exiting
    //****************************************************
 
    exit(0);
}
