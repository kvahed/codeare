/**************
 ** includes **
 **************/
#include <ctime>
#include <iostream>

# ifndef _timer_h_
#define _timer_h_

/******************
 ** declarations **
 ******************/
class MyTimer;


/*****************
 ** definitions **
 *****************/
// enums //
typedef enum {
        time_default,
        time_kernel,
        time_memcpy_HostToDevice,
        time_memcpy_DeviceToHost,
        time_check_config,
        time_init_malloc
} Measurement_type;

const char * const enum_names [] = { "time_default", "time_kernel", "time_memcpy_HostToDevice", "time_memcpy_DeviceToHost", "time_check_config", "time_init_malloc" };


class MyTimer {

private:
  static const int count = 6;
  bool indicator [count];
  clock_t mp_measurements [count];
  
public:
  double tic (const Measurement_type timer) {
    if (timer < count) {
      clock_t time = mp_measurements [timer];
      mp_measurements [timer] = clock ();
//      std::cout << "-----\"" << enum_names [timer] << "\"----- first / second : " << time << " / " << mp_measurements [timer] << ", value: " << ((double)(mp_measurements [timer] - time)) / CLOCKS_PER_SEC << std::endl;
      return (indicator [timer] ? indicator [timer] = false, .1 : ((double)(mp_measurements [timer] - time)) / CLOCKS_PER_SEC);
    } else {
      return -1;
    }
  }

};

#endif // ifndef _timer_h_
