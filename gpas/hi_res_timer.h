#ifndef _HI_RES_TIMER_H_
#define _HI_RES_TIMER_H_

#include <sys/time.h>

class HiResTimer
{
public:
    void start();
    void stop();
    unsigned long getElapsedTime();
private:
    timeval t1, t2;
};

#endif
